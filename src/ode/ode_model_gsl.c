#include "ode_model.h"

/* gsl_odeiv2_step_type *gsl_odeiv2_step_msbdf */
/* gsl_odeiv2_step_type *gsl_odeiv2_step_msadams */

typedef int (*dfdp) (double t, const double y[], double * dfdy, double dfdt[], void * params);
typedef int (*out_f) (double t, const double y[], double * yf, void * params);

/* also known as "ode_model"*/
struct solver_specific_ode_model {
  gsl_odeiv2_system *sys;
	yf *func;
  dfdp *jacp;
};

void *load_or_abort(void *lib, char *sym){
    fprintf(stderr,"[%s] loading «%s»",__func__,sym);
		void *s=dlsym(L,sym);
		assert(s);
		return s;
}

ode_model* ode_model_loadFromFile(const char *filename){
  void* L=dlopen(filename,RTLD_LAZY|RTLD_LOCAL);
  char *name=basename(filename);
  char *dot=strchr(name);
  int n=dot?dot-name:0;
  char *sym=alloca(n+12);
	ode_model *model = malloc(sizeof(ode_model));
  assert(system);
  memcpy(sym,name,n);
  if (L) {
    /* right hand side: f */
    model->sys=load_or_abort(L,strcpy(sym+n,"_sys"))
    /* df/dp */
    model->jacp=load_or_abort(L,strcpy(sym+n,"_jacp"));
		/* user defined functions */
		model->func=load_or_abort(L,strcpy(sym+n,"_out"));
  } else {
    fprintf(stderr, "Library %s could not be loaded.\n",filename);
    char *err = dlerror();
    if (err) fprintf(stderr, "%s",err);
    abort();
  }
	return model;
} 


ode_solver* ode_solver_alloc(ode_model *model){


}
