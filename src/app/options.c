#include "options.h"
#include "../ode/ode_model.h"
#include <assert.h>

void append(vec_d *a, double val){
  if (a->size==a->max){
    a->max+=4;
    a->value=realloc(a->value,sizeof(double)*(a->max));
  }
  a->value[a->size++]=val;
}

void clear(vec_d *a){
  a->size=0;
}

vec_d* vector_d(size_t n){
  vec_d *a=malloc(sizeof(vec_d));
  a->max=n;
  a->size=0;
  a->value=malloc(sizeof(double)*(a->max));
  return a;
}

void read_abs_tol(vec_d *a, char *abs_tol){
  assert(abs_tol);
  assert(a);
  clear(a);
  char *s=NULL,*p=abs_tol;
  double d;
  while (s!=p) {
    s=p;
    d=strtod(s,&p);
    if (s!=p) {
      append(a,d);
    }
  }
}

/* this is where the hard coded defaults are set. */
main_options* /* a struct with default values */
default_options(char *global_sample_filename_stem,/*for mcmc result files*/
		char *lib_name)/*Model library name*/{
  main_options *option=malloc(sizeof(main_options));
  option->initial_stepsize_rank_factor=1.0;
  option->output_file=global_sample_filename_stem;
  option->library_file=lib_name;
  option->target_acceptance=-0.25;
  option->initial_stepsize =-0.1;
  option->sample_size=-100;
  option->abs_tol=vector_d(10);
  append(option->abs_tol,ODE_SOLVER_ABS_ERR);
  option->rel_tol=ODE_SOLVER_REL_ERR;
  return option;
}

