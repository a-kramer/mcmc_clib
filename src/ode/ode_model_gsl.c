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

ode_ivp* ode_ivp_alloc(ode_model *M)
{
	assert(M);
	ode_ivp *solver_for_gsl=malloc(sizeof(ode_ivp));
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

void ode_ivp_free(ode_ivp *solver_for_gsl)
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

void ode_ivp_init(ode_ivp *s, const double t0, gsl_vector *y0, gsl_vector *p)
{
	gsl_vector_memcpy(s->params,p);
	s->sys.params=s->params->data;
	s->t0=t0;
	gsl_odeiv2_driver *drv=gsl_odeiv2_driver_alloc_y_new(s->odeModel->sys, gsl_odeiv2_step_msbdf, 1e-5, 1e-7, 1e-7);
	s->driver=drv;
}

/* solve A*X=B, all matrices */
int solve_linear_system(gsl_matrix *A, gsl_matrix *X, gsl_matrix *B){
	double *C;
	assert(A && X && B);
	size_t n=B->size1, m=B->size2;
	gsl_matrix_view x,b;
	gsl_matrix *T=gsl_matrix_alloc(n,n);

	int status=GSL_SUCCESS;
	if (A->size2 == X->size1){
		// C=alloca(sizeof(double)*n*m);
		// qr=gsl_matrix_view_array(C,n,m);
		status=gsl_linalg_QR_decomp_r(A,T);
		if (status!= GSL_SUCCESS){
			fprintf(stderr,"[%s] QR decomposition failed. %s\n",__func__,gsl_strerror(status));
		}
		for (i=0;i<m;i++){
			b=gsl_matrix_column(B,i);
			x=gsl_matrix_column(X,i);
			status|=gsl_linalg_QR_solve_r(A, T, &(b.vector), &(x.vector));
		}
	}
	gsl_matrix_free(T);
	return status;
}

static int approximate_sens(ode_model *m, double ti, gsl_vector *y_i, double tf, gsl_vector *y_f, gsl_matrix *yS){
	assert(yS);
	size_t ny=yS->size1;
	size_t np=yS->size2;
	assert(m->sys->dimesnion == ny);
	gsl_matrix *A=gsl_matrix_alloc(ny,ny);
	gsl_matrix *B=gsl_matrix_alloc(ny,np);
	double dfdt[ny];
  gsl_matrix_view Ayy=gsl_matrix_view_array(A,ny,ny);
  gsl_matrix_view Byp=gsl_matrix_view_array(A,ny,np);

  gsl_matrix *R;
  gsl_matrix *eJt;
  gsl_matrix *Jt=gsl_matrix_alloc(ny,ny);
  gsl_vector *tau=gsl_vector_alloc(ny);
  gsl_vector *x=gsl_vector_alloc(np);
  gsl_vector *r=gsl_vector_alloc(np);
	/**/
	gsl_vector_view eJt_diag;
	gsl_vector_view jacp_row;
	int status;

	int i;

	m->sys->jacobian(tf,y_f->data,A,dfdt,m->sys->params);
	m->jacp(tf,y_f->data,B,dfdt,m->sys->params);

	//gsl_matrix_transpose(a->jacobian_y); // now jacobian_y(i,j)=df[i]/dy[j];
	gsl_matrix_memcpy(Jt,A);
	gsl_matrix_scale(Jt,tf-ti);
	gsl_matrix *L=gsl_matrix_alloc(ny,np); /* jac_y\jac_p in matlab syntax */
  solve_linear_system(A,L,B);

	gsl_matrix *R=gsl_matrix_alloc(n,m);
	gsl_matrix_memcpy(R,L);
	gsl_matrix_add(R,yS);
	status=gsl_linalg_exponential_ss(Jt,eJt,GSL_PREC_SINGLE);
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

static int transform(ode_vector *v, struct tf *tfm)
{
	if (tfm && a){
		switch (tfm->type){
		case tf_matrix_vector:
			gsl_vector_memcpy(tfm->result,tfm->b);
			gsl_blas_dgemv(CblasNoTrans, 1.0, tfm->A, v, 1.0, tfm->result); /* A v + b */
			gsl_vector_memcpy(v,tfm->result);
			break;
		case tf_vector_vector:
			gsl_vector_mul(v,tfm->a);
			gsl_vector_add(v,tfm->b);
			break;
		case tf_vector_scalar:
			gsl_vector_mul(v,tfm->a);
			gsl_vector_add_constant(v,tfm->s_b)
			break;
		case tf_scalar_scalar:
			gsl_vector_scale(v,tfm->s_a);
			gsl_vector_add_constant(v,tfm->s_b);
			break;
		}
	}
	return GSL_SUCCESS;
}

/* this functions solves the initial value problem and model in s,
	 using the driver struct in s; the result is returned in y and fy; the
	 solver contains pre-allocated memory to store the model state
	 during integration */
int ode_ivp_advance(ode_ivp *s, gsl_vector *t, gsl_vector **y, gsl_vector **fy, gsl_matrix **yS, gsl_matrix **fyS, struct scheduled_event **e)
{
	int j;
	double ti=s->t0; /* initial */
	double tf;       /* final */
	int status;
	struct scheduled_evemt *ej;
	gsl_odeiv2_driver_reset(s->driver);
	for (j=0; j<t->size; j++){
		tf=gsl_vector_get(t,j);
		/* loop over events prior to t[j] */
		while (e && e[j] && e[j]->t < tf){
			status=gsl_odeiv2_driver_apply(s->driver, &ti, e[j]->t, s->y->data);
			transform(s->y,e[j]->state);
			transform(s->params,e[j]->params);
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

