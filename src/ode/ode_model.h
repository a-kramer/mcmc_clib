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
	void *driver;
	ode_model *odeModel;
	gsl_vector *y;
	gsl_matrix *yS;
	gsl_vector *fy;
	gsl_matrix *fyS;
	gsl_matrix *jac;
	gsl_matrix *jacp;
	gsl_vector *params;
} ode_solver;

/* Loads shared library with user defined functions and ode model data */
ode_model* ode_model_loadFromFile(const char *filename);
/* Creates a new ode_solver for the given model */
ode_solver* ode_solver_alloc(ode_model* model);
/* Initialises the ode solver. y0, yS0 and p can be NULL. */
void ode_solver_init(ode_solver* solver, const double t0, gsl_vector *y0, gsl_vector *p);
/* Re-initialises the ode solver. y0, yS0 and p can be NULL in which case the default values are used.*/
void ode_solver_reinit(ode_solver* solver, const double t0,  double* y0, int lenY, const double* p, int lenP );
int ode_solver_solve(ode_solver* solver, double *t, double tf);
/* Returns sensitivities for the current solution at t. This function must be called only after ode_solver_solve. */
void ode_solver_free(ode_solver* solver);
#endif
