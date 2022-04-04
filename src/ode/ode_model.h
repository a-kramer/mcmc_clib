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

enum tf_type = {tf_matrix_vector, tf_vector_vector, tf_vector_scalar, tf_scalar_scalar};

typedef struct sens_approx sensitivity_approximation;

typedef struct {
	void *driver;
	ode_model *odeModel;
	double t0;
	gsl_vector *y;
	gsl_matrix *yS;
	gsl_vector *fy;
	gsl_matrix *fyS;
	gsl_matrix *jac;
	gsl_matrix *jacp;
	gsl_vector *params;
	sensitivity_approximation *sapprox;
} ode_solver;

struct tf {
	tf_type type;
	union {
		gsl_matrix *A;
		gsl_vector *a;
		double s_a;
	};
	union {
		gsl_vector *b;
		double s_b;
	};
	gsl_vector *result;
};

struct scheduled_event {
	double t;
	struct tf *state;
	struct tf *params;
	struct scheduled_event *next;
};

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
