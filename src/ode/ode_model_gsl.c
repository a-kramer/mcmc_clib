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
	gsl_odeiv2_driver *drv=gsl_odeiv2_driver_alloc_y_new(s->odeModel->sys, gsl_odeiv2_step_msbdf, 1e-5, 1e-7, 1e-7);
	s->driver=drv;
}

//void ode_solver_reinit(ode_solver* solver, const double t0,  double* y0, int lenY, const double* p, int lenP );
void ode_solver_reinit(ode_solver *s, const double t0, gsl_vector *y0, gsl_vector *p)
{
	gsl_odeiv2_driver_reset(s->driver);
}

int ode_solver_solve(ode_solver* s, double *t, const double tf)
{
	return gsl_odeiv2_driver_apply(s->driver, t, tf, s->y->data);
}
