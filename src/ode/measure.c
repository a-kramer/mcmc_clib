#include <gsl/gsl_vector.h>
#include "../ode/ode_solver.h"

int measure(char *lib_name, gsl_vector *t, realtype *solver_param){

  /* load model */
  ode_model* odeModel = ode_model_loadFromFile(lib_name);  /* alloc */
  if (odeModel == 0) {
    fprintf(stderr, "Library %s could not be loaded.\n",lib_name);
    exit(1);
  }
  
  
  /* define local variables for parameters and inital conditions */
  int N = ode_model_getN(odeModel);
  int P = ode_model_getP(odeModel);
  int F = ode_model_getF(odeModel);
  const char** v_names = ode_model_get_var_names(odeModel);
  const char** p_names = ode_model_get_param_names(odeModel);
  const char** f_names = ode_model_get_func_names(odeModel);
  
  double y[N];
  ode_model_get_initial_conditions(odeModel, y, N);
  
  double p[P];
  ode_model_get_default_params(odeModel, p, P);

  /* create solver */
  ode_solver* solver = ode_solver_alloc(odeModel); /* alloc */

  if (solver == 0) {
    fprintf(stderr, "Solver %s could not be created.\n",lib_name);
    ode_model_free(odeModel);
    if(outputF != stdout)
      fclose(outputF);
    exit(1);
  }

  /* init solver */
  double t_old = solver_param[2];
  ode_solver_init(solver, t_old, y, N, p, P);
  ode_solver_setErrTol(solver, solver_param[1], &solver_param[0], 1);
  if (ode_model_has_sens(odeModel)) {
    ode_solver_init_sens(solver, NULL, 0, 0);
  }
  
  double tout;
  for (i=0; i < t->size; i++){

    int CVerror =  ode_solver_solve(solver, gsl_vector_get(t,i), y, &tout);
    if (CVerror) {
      fprintf(stderr, "ODE solver failed.\n");
      goto cleanup;
    }
    
    PrintOutputVars(outputF, tout, y, N, 1);
    
    if (ode_model_has_funcs(odeModel)) {
      double fy[F];
      ode_solver_get_func(solver, tout, y, fy);
      PrintOutputVars(outputF, tout, fy, F, 0);
    }
    
    if(ode_model_has_sens(odeModel)){
      double yS[P*N];
      ode_solver_get_sens(solver, tout, yS);
      PrintOutputSens(outputF, tout, yS, P, N, 0);
      
      if (ode_model_has_funcs_sens(odeModel)) {
	double fyS[P*F];
	ode_solver_get_func_sens(solver, tout, y, yS, fyS);
	PrintOutputSens(outputF, tout, fyS, P, F, 0);
      }
    }
  }
  
  
  
}
