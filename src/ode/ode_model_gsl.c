#include "ode_model.h"

/* gsl_odeiv2_step_type *gsl_odeiv2_step_msbdf */
/* gsl_odeiv2_step_type *gsl_odeiv2_step_msadams */

typedef int (*dfdp) (double t, const double y[], double * dfdy, double dfdt[], void * params);
typedef int (*out_f) (double t, const double y[], double * yf, void * params);

/* also known as "ode_model" */
struct solver_specific_ode_model {
	gsl_odeiv2_system sys;
	dfdp *jacp;
	out_f *func;
};

/* also known as sensitivity_approximation */
struct sens_approx {
  gsl_matrix *temp_ny_ny;
  gsl_matrix *temp_ny_np;
  double t0;
  gsl_matrix *R;
  gsl_matrix *eJt;
  gsl_matrix *Jt;
  gsl_vector *tau;
  gsl_vector *x;
  gsl_vector *r;
};

void *load_or_abort(void *lib, char *sym)
{
	fprintf(stderr,"[%s] loading «%s»",__func__,sym);
	void *s=dlsym(L,sym);
	assert(s);
	return s;
}

/* loads functions from shared library, all functions are necessary, returned struct can be freed by free() */
ode_model* ode_model_loadFromFile(const char *filename)
{
	void* L=dlopen(filename,RTLD_LAZY|RTLD_LOCAL);
	char *name=basename(filename);
	char *dot=strchr(name,'.');
	size_t n=dot?dot-name:0;
	int d;
	char *sym=alloca(n+12);
	// char *suffix[]={"_vf", "_jac", "_jacp", "_out"};
	// int n=sizeof(suffix)/sizeof(char*);
	assert(sym);
	ode_model *model = malloc(sizeof(struct solver_specific_ode_model));
	int i;
	memcpy(sym,name,n);
	if (L) {
		/* ode right-hand-side function f (a vector field) */
		model->sys=malloc(sizeof(gsl_odeiv2_system));
		model->sys.function=load_or_abort(L,strcpy(sym+n,"_vf"));
		/* try evaluation, vf must accept NULL pointers and return y-size */
		d=model->sys.function(0.0,NULL,NULL,NULL);
		assert(d>0);
		model->sys.dimension=d;
		/* df/dy */
		model->sys.jacobian=load_or_abort(L,strcpy(sym+n,"_jac"));
		assert(model->sys.jacobian(0.0,NULL,NULL,NULL,NULL)==d*d);
		/* df/dp */
		model->jacp=load_or_abort(L,strcpy(sym+n,"_jacp"));
		/* user defined output functions */
		model->func=load_or_abort(L,strcpy(sym+n,"_out"));
	} else {
		fprintf(stderr, "Library %s could not be loaded.\n",filename);
		char *err = dlerror();
		if (err) fprintf(stderr, "%s",err);
		abort();
	}
	return model;
}

ode_solver* ode_solver_alloc(ode_model *M)
{
	assert(M);
	ode_solver *solver_for_gsl=malloc(sizeof(ode_solver));
	int ny=M->sys.dimension;
	int nfy=M->func(0.0,NULL,NULL,NULL);
	int npar=M->dfdp(0.0,NULL,NULL,NULL)/ny;
	solver_for_gsl->odeModel=M;
	solver_for_gsl->y=gsl_vector_alloc(ny);
	solver_for_gsl->fy=gsl_vector_alloc(nfy);
	solver_for_gsl->params=gsl_vector_alloc(npar);
	solver_for_gsl->jac=gsl_matrix_alloc(ny,ny);
	solver_for_gsl->jacp=gsl_matrix_alloc(npar,npar);
	solver_for_gsl->yS=gsl_matrix_alloc(ny,npar);
	solver_for_gsl->fyS=gsl_matrix_alloc(nfy,npar);
	return solver_for_gsl;
}

void ode_solver_free(ode_solver *solver_for_gsl)
{
	gsl_vector_free(solver_for_gsl->y);
	gsl_vector_free(solver_for_gsl->fy);
	gsl_vector_free(solver_for_gsl->params);
	gsl_matrix_free(solver_for_gsl->jac);
	gsl_matrix_free(solver_for_gsl->jacp);
	gsl_matrix_free(solver_for_gsl->yS);
	gsl_matrix_free(solver_for_gsl->fyS);
	gsl_odeiv2_driver_free(solver_for_gsl->driver);
	free(solver_for_gsl->odeModel);
	free(solver_for_gsl);
}

void ode_solver_init(ode_solver *s, const double t0, gsl_vector *y0, gsl_vector *p)
{
	gsl_vector_memcpy(s->params,p);
	s->sys.params=s->params->data;
	s->t0=t0;
	gsl_odeiv2_driver *drv=gsl_odeiv2_driver_alloc_y_new(s->odeModel->sys, gsl_odeiv2_step_msbdf, 1e-5, 1e-7, 1e-7);
	s->driver=drv;
}


static int approximate_sens(ode_solver *solver, double tout, gsl_vector *y, gsl_vector *fy, gsl_matrix *yS, gsl_matrix *fyS, sensitivity_approximation *a){
	gsl_vector_view eJt_diag;
	gsl_vector_view jacp_row;
	int status;
	assert(yS);
	int P=yS->size1;
	int i;
	assert(a && a->jacobian_y && a->jacobian_p && a->tau && a->x && a->r);
	assert(y && fy && yS && fyS);

	/* CVODE stores matrices column wise, but gsl stores them row wise */
	/* so copying values from CVODE to gsl transposes a matrix */
	ode_solver_get_jac(solver,tout,y->data,fy->data,a->jacobian_y->data);
	ode_solver_get_jacp(solver,tout,y->data,fy->data,a->jacobian_p->data);

	//gsl_matrix_transpose(a->jacobian_y); // now jacobian_y(i,j)=df[i]/dy[j];
	gsl_matrix_memcpy(a->Jt,a->jacobian_y);
	gsl_matrix_scale(a->Jt,tout-(a->t0));				// this is now Jacobian*t
	status=gsl_linalg_QR_decomp(a->jacobian_y,a->tau);
	if (status!= GSL_SUCCESS){
		fprintf(stderr,"[%s] QR decomposition failed. %s\n",__func__,gsl_strerror(status));
	}
	fflush(stdout);
	for (i=0;i<P;i++){
		jacp_row=gsl_matrix_row(a->jacobian_p,i);
		// solve in place; jac_p will contain the solution: jac_y\jac_p
		status=gsl_linalg_QR_svx(a->jacobian_y, a->tau, &jacp_row.vector);
		if (status!=GSL_SUCCESS) {
			fprintf(stderr,"[%s] QR solution of linear equations failed: %s. Using minimum Norm solution.\n",__func__,gsl_strerror(status));
			if (gsl_linalg_QR_lssolve(a->jacobian_y, a->tau, &jacp_row.vector, a->x, a->r)!=GSL_SUCCESS) {
				fprintf(stderr,"[%s] QR lssolve also failed. exiting.\n",__func__);
				abort();
			} else {
				gsl_vector_memcpy(&jacp_row.vector,a->x);
			}
		}
	}
	gsl_matrix *L=a->jacobian_p; /* jac_y\jac_p in matlab syntax */
	//gsl_matrix *R=gsl_matrix_alloc(s1,s2);
	gsl_matrix_memcpy(a->R,L);
	gsl_matrix_add(a->R,yS);
	status=gsl_linalg_exponential_ss(a->Jt,a->eJt,GSL_PREC_SINGLE);
	if (status!=GSL_SUCCESS){
		// this is not yet considered stable by GSL :/
		fprintf(stderr,"[%s] matrix exponential failed. %s\n",__func__,gsl_strerror(status));
		return GSL_FAILURE;
	}
	eJt_diag=gsl_matrix_diagonal(a->eJt);
	gsl_matrix_memcpy(yS,L);
	gsl_vector_add_constant(&eJt_diag.vector,-1.0);
	if (gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0,a->R,a->eJt,-1.0,yS)!=GSL_SUCCESS){
		fprintf(stderr,"[%s] general double precision	Matrix·Matrix failed.\n",__func__);
		return GSL_FAILURE;
	}
	return GSL_SUCCESS;
}

static int transform(ode_vector *a, struct tf *tfm)
{
	if (tfm && a){
		switch (tfm->type){
		case tf_matrix_vector:
			gsl_vector_memcpy(tfm->result,tfm->b);
			gsl_blas_dgemv(CblasNoTrans, 1.0, tfm->A, a, 1.0, tfm->result);
			gsl_vector_memcpy(a,tfm->result);
			break;
		case tf_vector_vector:
			gsl_vector_mul(a,tfm->a);
			gsl_vector_add(a,tfm->b);
			break;
		case tf_vector_scalar:
			gsl_vector_mul(a,tfm->a);
			gsl_vector_add_constant(a,tfm->s_b)
			break;
		case tf_scalar_scalar:
			gsl_vector_scale(a,tfm->s_a);
			gsl_vector_add_constant(a,tfm->s_b);
			break;
		}
	}
	return GSL_SUCCESS;
}

/* this functions solves the initial value problem and model in s,
	 using the driver struct in s; the result is returned in y and fy; the
	 solver contains pre-allocated memory to store the model state
	 during integration */
int ode_solver_apply(ode_solver *s, gsl_vector *t, gsl_vector **y, gsl_vector **fy, gsl_matrix **yS, gsl_matrix **fyS, struct scheduled_event **e)
{
	int j;
	double ti=s->t0; /* initial */
	double tf;       /* final */
	int status;
	struct scheduled_evemt *ej;
	gsl_odeiv2_driver_reset(s->driver);
	for (j=0; j<t->size; j++){
		ej=e[j];
		tf=gsl_vector_get(t,j);
		/* loop over events prior to tf */
		while (ej && ej->t < tf){
			status=gsl_odeiv2_driver_apply(s->driver, &ti, ej->t, s->y->data);
			transform(s->y,ej->state);
			transform(s->params,ej->state);
			assert(status==GSL_SUCCESS);
		}
		/* now advance to the next save worthy time point */
		status=gsl_odeiv2_driver_apply(s->driver, &ti, tf, s->y->data);
		switch (status){
		case GSL_SUCCESS:
			gsl_vector_memcpy(y[j],s->y);
			if (fy && fy[j]) s->odeModel->func(tf,y[j]->data,fy[j]->data,s->sys->params);
			if (yS) {
				approximate_sens(s, ti, a);
			}
			break;
		default:
			fprintf(stderr,"[%s] integration error: %i\n",__func__,status);
		}
	}
	return j;
}

int ode_solver_process_sens(
	ode_solver *solver,
	double tout,
	gsl_vector *y, gsl_vector* fy,
	gsl_matrix *yS, gsl_matrix *fyS,
	sensitivity_approximation *a)
{
	ode_model *model;
	model=solver->odeModel;

	if (ode_model_has_sens(model)){
		ode_solver_get_sens(solver, tout, yS->data);
	} else {
		approximate_sens(solver,tout,y,fy,yS,fyS,a);
	}
	if (ode_model_has_funcs_sens(model)){
		ode_solver_get_func_sens(solver, tout, y->data, yS->data, fyS->data);
	}
	return GSL_SUCCESS;
}

