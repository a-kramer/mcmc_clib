/*
 *  ode_model.h
 *  odeSolve
 *
 *  Created by Vassilios Stathopoulos on 28/10/2011.
 *  Copyright 2011 Computing Science. All rights reserved.
 *
 */
#ifndef ODE_MODEL_H
#define ODE_MODEL_H
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct ivp_specific_ode_model ode_model;

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
	gsl_matrix *fyJ;
	gsl_matrix *fyP;
	gsl_vector *params;
	sensitivity_approximation *sapprox;
} ode_ivp;

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
/* Creates a new ode_ivp for the given model */
ode_ivp* ode_ivp_alloc(ode_model* model);
/* Initialises the ode ivp. y0, yS0 and p can be NULL. */
void ode_ivp_init(ode_ivp* ivp, const double t0, gsl_vector *y0, gsl_vector *p);
/* Re-initialises the ode ivp. y0, yS0 and p can be NULL in which case the default values are used.*/
int ode_ivp_advance(ode_ivp *s, gsl_vector *t, gsl_vector **y, gsl_vector **fy, gsl_matrix **yS, gsl_matrix **fyS, struct scheduled_event **e);
/* Returns sensitivities for the current solution at t. This function must be called only after ode_ivp_solve. */
void ode_ivp_free(ode_ivp* ivp);
#endif
