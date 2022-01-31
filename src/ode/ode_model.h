/*
 *  ode_model.h
 *  odeSolve
 *
 *  Created by Vassilios Stathopoulos on 28/10/2011.
 *  Copyright 2011 Computing Science. All rights reserved.
 *
 */
#ifndef _ODE_MODEL_H
#define _ODE_MODEL_H
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct solver_specific_ode_model ode_model;

typedef struct{
  void* solver_mem;
  ode_model*  odeModel;
  gsl_vector* y;
  gsl_matrix* yS;
  gsl_vector* fy;
	gsl_matrix* fyS;
  gsl_matrix* jac;
  gsl_matrix* jacp;
  gsl_vector* params;
} ode_solver;

  
/* Loads shared library with user defined functions and ode model data */
ode_model* ode_model_loadFromFile(const char *filename);

/* returns a constant pointer to a variable names array */
const char** ode_model_get_var_names(const ode_model* model);
	
/* returns a constant pointer to a parameter names array */
const char** ode_model_get_param_names(const ode_model* model);
	
/* returns a constant pointer to a function names array */
const char** ode_model_get_func_names(const ode_model* model);
  
/* Frees resources used by ode_model */
void ode_model_free(ode_model* model);

/* Creates a new ode_solver for the given model */
ode_solver* ode_solver_alloc(ode_model* model);
	
/* Initialises the ode solver. y0, yS0 and p can be NULL. */
void ode_solver_init(ode_solver* solver, const double t0, double* y0, int lenY, double* p, int lenP );
void ode_solver_init_sens(ode_solver* solver,  double* yS0, int lenP, int lenY);

/* Sets error tollerances, NOTE this function should be called only after an ode_solver_init or ode_solver_reinit. */
void ode_solver_setErrTol(ode_solver* solver, const double rel_tol, double* abs_tol, const int abs_tol_len);
	
/* Re-initialises the ode solver. y0, yS0 and p can be NULL in which case the default values are used.*/
void ode_solver_reinit(ode_solver* solver, const double t0,  double* y0, int lenY, const double* p, int lenP );
void ode_solver_reinit_sens(ode_solver* solver, double* yS0, int lenP, int lenY);
	
/* Solves the ode system until time t and returns solution in y.*/	
int ode_solver_solve(ode_solver* solver, const double t, double* tout);
/* Returns sensitivities for the current solution at t. This function must be called only after ode_solver_solve. */	
void ode_solver_get_sens(ode_solver* solver, double t, double* yS);
/* Returns functions of the current solution. */	
void ode_solver_get_func(ode_solver* solver, const double t, double* y, double* fy);
/* Returns sensitivities of functions of the current solution. */	
void ode_solver_get_func_sens(ode_solver* solver, const double t, double* y, double* yS, double* fyS);
void ode_solver_get_jac(ode_solver* solver, const double t,  double* y,  double* fy, double *jac);
void ode_solver_get_jacp(ode_solver* solver, const double t,  double* y,  double* fy, double *jacp);
/* Prints to out solver statistics. */
/* void ode_solver_print_stats(const ode_solver* solver, FILE* outF); */
/* Frees resources used by the ode_solver */
void ode_solver_free(ode_solver* solver);
#endif
