#ifndef MODEL_PARAMETERS_H
#define MODEL_PARAMETERS_H
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <stdlib.h>
#include <stdio.h>
#include "../ode/ode_model.h"
#include "../app/event.h"
// prior types
#include "ptype.h"
// normalisation methods
#define DATA_IS_ABSOLUTE 0
#define DATA_NORMALISED_BY_REFERENCE 1
#define DATA_NORMALISED_BY_TIMEPOINT 2
#define DATA_NORMALISED_BY_STATE_VAR 3
#define DATA_NORMALISED_INDIVIDUALLY 4
// for individual normalisation:
#define NEEDS_NORMALISATION(E) ((E->NormaliseByExperiment>=0) || (E->NormaliseByTimePoint)>=0)

/* This structure contains some pre-allocated space for normalisation
 * purposes.
 */
struct normalization{
	gsl_vector **fy; // normalising fy; simulation result
	gsl_matrix **fyS; // normalising fyS; simulation result
	gsl_vector **data; // normalising data vector	(from a lab)
	gsl_vector **stdv; // and its standard deviation;
};


/* Data can be stored in one of two places:
 * 1. in the model parameters structure (in the Data matrix).
 * 2. in the experiment structure, individually for each experiment.
 *
 * If method 1 is used, then the experiments have a view of the part
 * of [Data] that belongs to them. Otherwise data_block is just a normal matrix.
 *
 * If method 2 is used, then the main Data matrix is NULL, and
 * data_block is actually the primary owner of its contents.
 */
typedef struct {
	double t0;
	gsl_vector *t;
	struct normalization *normalise;
	double pdf_lognorm; // normalisation constant of the probability density
	int lflag;								// TRUE if this data should influence LogLikelihood explicitely; FALSE if it is used in some implicit way (e.g. normalisation)
	gsl_vector **data;				// data at time[j] is accessed as experiment[i]->data[j]
	gsl_vector **sd_data;
	gsl_vector **y; // y[t](i) vector of size (t->size);
	gsl_vector **fy; // measurement model output functions
	gsl_vector *init_y;
	gsl_vector *input_u;
	gsl_matrix *yS0;
	gsl_matrix **yS;
	gsl_matrix **fyS;
	gsl_matrix **oS; // observational sensitivities
	int NormaliseByExperiment;
	int NormaliseByTimePoint;
	gsl_vector_int *NormaliseByOutput;
	event_list_t *event_list;
	event_row_t **single; // array of linked lists
	before_measurement **before_t; // array of pointers, converted from linked list
	int index[3];
} experiment;

typedef struct {
	union{
		gsl_vector *mu;
		gsl_vector *alpha;
	};
	gsl_vector **tmp; // temporary storage for prior calculations
	int n; // number of temporary vectors above
	union {
		gsl_matrix *Sigma_LU;	// LU factors of the Sigma matrix
		gsl_matrix *inv_cov;	 // Same matrix already given in inverted form
		gsl_vector *sigma;	// If Sigma is diagonal, this is the diagonal vector.
		gsl_vector *beta;	// this is a parameter of generalised Gaussians
	};
	gsl_permutation *p; // used whenever LU decomposition is used
	int signum;				 // also belongs to LU decomposition
	int type;					 // for later use if maybe generalised prior is available, but also, diagonal/non-diagonal stuff;
} prior_t;

typedef struct {
	double t0;
	problem_size *size;
	double pdf_lognorm;
	gsl_vector *p; // memory for the ODE's parameters (they are
		 // exponential) and input parameters, appended: p=[exp(x),u];
	experiment **E; // an experiment (data, initial conditions, inputs, simulation results)
	int normalisation_type;
	gsl_vector *tmpF; // temporay storage of size F
	gsl_matrix *tmpDF; // temporary storage of size D×F
	prior_t *prior;
	sensitivity_approximation **S_approx;
	ode_solver **solver; // contains: cvode_mem; *odeModel; N_Vector y; N_Vector *yS; params; per experiment
} bayesian_parameters;

int bayesian_parameters_alloc(bayesian_parameters*);
int bayesian_parameters_free(bayesian_parameters*);
int bayesian_parameters_init(bayesian_parameters*);
#endif
