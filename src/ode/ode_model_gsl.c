#include "ode_model.h"
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_linalg.h>
#include <dlfcn.h>
#include <assert.h>
#include <string.h>
#include <libgen.h>
/* gsl_odeiv2_step_type *gsl_odeiv2_step_msbdf */
/* gsl_odeiv2_step_type *gsl_odeiv2_step_msadams */

typedef int (*dfdp) (double t, const double y[], double * dfdy, double dfdt[], void * params);
typedef int (*out_f) (double t, const double y[], double * F, void * params);
typedef int (*out_f_dy) (double t, const double y[], double * dFdy, void * params);
typedef int (*out_f_dp) (double t, const double y[], double * dFdp, void * params);
typedef int (*f_default_parameters) (double t, void * params);
typedef int (*f_initial_condition) (double t, double *y, void * params);

/* also known as "ode_model" */
struct solver_specific_ode_model {
	gsl_odeiv2_system sys;
	f_default_parameters param_default;
	f_initial_condition init;
	dfdp jacp;
	out_f func;
	out_f_dy func_jac;
	out_f_dp func_jacp;
};

void *load_or_abort(void *lib, char *sym)
{
	fprintf(stderr,"[%s] loading «%s»",__func__,sym);
	void *s=dlsym(lib,sym);
	assert(s);
	return s;
}

void *load_optional(void *lib, char *sym)
{
	fprintf(stderr,"[%s] loading «%s»",__func__,sym);
	void *s=dlsym(lib,sym);
	return s;
}


/* loads functions from shared library, all functions are necessary, returned struct can be freed by free() */
ode_model* ode_model_load_from_file(char *filename)
{
	void* L=dlopen(filename,RTLD_LAZY|RTLD_LOCAL);
	char *name=basename(filename);
	char *dot=strchr(name,'.');
	size_t n=dot?dot-name:0;
	int d;
	char *sym=alloca(n+32);
	// char *suffix[]={"_vf", "_jac", "_jacp", "_out"};
	// int n=sizeof(suffix)/sizeof(char*);
	assert(sym);
	ode_model *model = malloc(sizeof(struct solver_specific_ode_model));
	memcpy(sym,name,n);
	if (L) {
		/* ode right-hand-side function f (a vector field) */
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
		model->param_default=load_optional(L,strcpy(sym+n,"_default"));
		model->init=load_optional(L,strcpy(sym+n,"_init"));
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
	int npar=M->jacp(0.0,NULL,NULL,NULL,NULL)/ny;
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

/* initializes the ode_ivp struct with values taken from t0, y0, and
	 p. The memory of y0 and p may be freed afterwards. if y0 is NULL
	 this function will use the model's defaults (_init function). If
	 the parameters p are NULL the models _default function will be used
	 to set default parameters.  Thus, `ode_ivp_init(s,0.0,NULL,NULL)`
	 is a valid call. In cases where the _default or _init function is
	 missing from the model, all values are set to 1.0 instead.
*/
void ode_ivp_init(ode_ivp *s, const double t0, gsl_vector *y0, gsl_vector *p)
{
	s->t0=t0;
	if (p){
		gsl_vector_memcpy(s->params,p);
	} else if (s->odeModel->param_default != NULL){
		s->odeModel->param_default(t0,s->params);
	} else {
		gsl_vector_set_all(s->params,1.0);
	}
	s->odeModel->sys.params=s->params->data; /* gsl_odeiv2 solvers will use this pointer */
	if (y0) {
		gsl_vector_memcpy(s->y,y0);
	} else if (s->odeModel->init != NULL){
		s->odeModel->init(t0,s->y->data,s->params->data);
	} else {
		gsl_vector_set_all(s->y,1.0);
	}
	gsl_matrix_set_zero(s->yS);
	gsl_matrix_set_zero(s->fyS);
	gsl_odeiv2_driver *drv=gsl_odeiv2_driver_alloc_y_new(&(s->odeModel->sys), gsl_odeiv2_step_msbdf, 1e-5, 1e-7, 1e-7);
	s->driver=drv;
}

/* solve A*X=B, all matrices */
int solve_linear_system(gsl_matrix *A, gsl_matrix *X, gsl_matrix *B){
	assert(A && X && B);
	int i;
	size_t n=B->size1, m=B->size2;
	gsl_vector_view x,b;
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

	/* this is according to this approximation:
	yS0=yS(t0);
	yS(t) = exp(A*(t-t0))*(yS0 + A\B) - A\B
	eJt = exp(A*(t-t0))
	C = A\B
	R = yS0 + C
	yS = eJt*R - C
	*/
static int approximate_sens(ode_model *m, double ti, gsl_vector *y_i, double tf, gsl_vector *y_f, gsl_matrix *yS){
	assert(yS);
	size_t ny=yS->size1;
	size_t np=yS->size2;
	assert(m->sys.dimension == ny);
	gsl_matrix *A=gsl_matrix_alloc(ny,ny);
	gsl_matrix *B=gsl_matrix_alloc(ny,np);
	double dfdt[ny];
	gsl_matrix *C=gsl_matrix_alloc(ny,np);
	gsl_matrix *R=gsl_matrix_alloc(ny,np);
	gsl_matrix *eJt=gsl_matrix_alloc(ny,ny);
	gsl_matrix *Jt=gsl_matrix_alloc(ny,ny);
	int status;
	m->sys.jacobian(tf,y_f->data,A->data,dfdt,m->sys.params); /* A is df/dy (the jacobian) */
	m->jacp(tf,y_f->data,B->data,dfdt,m->sys.params);          /* B is df/dp (the parameter jacobian) */
	gsl_matrix_memcpy(Jt,A);
	gsl_matrix_scale(Jt,tf-ti);
	solve_linear_system(A,C,B); /* solve A*C=B, or equivalently C = A\B */
	gsl_matrix_memcpy(R,C);
	gsl_matrix_add(R,yS);
	if ((status=gsl_linalg_exponential_ss(Jt,eJt,GSL_PREC_SINGLE))!=GSL_SUCCESS){
		fprintf(stderr,"[%s] matrix exponential failed. %s\n",__func__,gsl_strerror(status));
	}
	gsl_matrix_memcpy(yS,C);
	if ((status=gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,1.0,R,eJt,-1.0,yS))!=GSL_SUCCESS){
		fprintf(stderr,"[%s] Matrix product (dgemm) failed. %s\n",__func__,gsl_strerror(status));
	}
	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_matrix_free(Jt);
	gsl_matrix_free(eJt);
	gsl_matrix_free(R);
	gsl_matrix_free(C);
	return status;
}

static int transform(gsl_vector *v, struct tf *tfm)
{
	if (tfm){
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
			gsl_vector_add_constant(v,tfm->s_b);
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
int ode_ivp_solve(ode_ivp *s, gsl_vector *t, gsl_vector **y, gsl_vector **fy, gsl_matrix **yS, gsl_matrix **fyS, struct scheduled_event **e)
{
	assert(s);
	int j;
	double ti=s->t0; /* initial */
	double tf;       /* final */
	int status;
	double *dfdt=alloca(s->y->size);
	assert(s->t0==gsl_vector_get(t,0));
	gsl_vector_memcpy(y[0],s->y);
	gsl_matrix_memcpy(yS[0],s->yS);
	gsl_odeiv2_driver_reset(s->driver);
	for (j=1; j<t->size; j++){
		tf=gsl_vector_get(t,j);
		/* loop over events prior to t[j] */
		while (e && e[j] && e[j]->t < tf){
			status=gsl_odeiv2_driver_apply(s->driver, &ti, e[j]->t, s->y->data);
			transform(s->y,e[j]->state);
			transform(s->params,e[j]->params);
			gsl_odeiv2_driver_reset(s->driver);
			assert(status==GSL_SUCCESS);
		}
		/* now advance to the next save worthy time point */
		status=gsl_odeiv2_driver_apply(s->driver, &ti, tf, s->y->data);
		switch (status){
		case GSL_SUCCESS:
			gsl_vector_memcpy(y[j],s->y);
			s->odeModel->sys.jacobian(tf,s->y->data,s->jac->data,dfdt,s->odeModel->sys.params);
			s->odeModel->jacp(tf,s->y->data,s->jacp->data,dfdt,s->odeModel->sys.params);
			if (fy && fy[j]) s->odeModel->func(tf,y[j]->data,fy[j]->data,s->odeModel->sys.params);
			if (yS) approximate_sens(s->odeModel, ti, y[j-1], tf, s->y, s->yS);
			if (yS && fyS && s->odeModel->func_jac && s->odeModel->func_jacp) {
				s->odeModel->func_jac(tf,y[j]->data,s->fyJ->data,s->odeModel->sys.params);
				s->odeModel->func_jacp(tf,y[j]->data,s->fyP->data,s->odeModel->sys.params);
			}
			break;
		default:
			fprintf(stderr,"[%s] integration error: %i\n",__func__,status);
		}
	}
	return j;
}

