#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdlib.h>
#include <stdio.h>
#include "../ode/ode_model.h"

typedef struct {
  int D; // number of sampling variables; unknown ODE-Model parameters
  int P; // number of total parameters of the ODE-Model (including dependent and known parameters)
  int N; // number of state varibales
  int F; // number of output functions
  int T; // number of measurement time-points for NORMALISATION_BY_REFERENCE
  int U; // number of input parameters
  int C; // number of experimental conditions, excluding control (reference measurement)
  int R; //number of reference measurements, either 1 or C
} problem_size;

typedef struct {
  double t0;
  gsl_vector *t;
  gsl_matrix *Data; // this is a block of size T*F for each experiment
  gsl_matrix *sdData; // standard deviation of the datapoints
  gsl_vector_view *data_row; // for convenience, we define some matrix views
  gsl_vector **data;          // now we can view any of the data rows as a vector
  gsl_vector_view *sd_data_row; // same for the standard deviations.
  gsl_vector **sd_data;
  gsl_vector **y; // y[c*T+j](i)
  gsl_vector **fy; // measurement model output functions
  gsl_vector *init_y;
  gsl_vector *input_u;
  gsl_matrix *yS0;
  gsl_matrix **yS;
  gsl_matrix **fyS;
  gsl_matrix **oS; // observational sensitivities
} experiment;

typedef struct {
  int D;
  double beta; // inverse temperature for annealing or tempering; logPosterior= beta*logLikelihood+logPrior
  gsl_vector *p; // memory for the ODE's parameters (they are
		 // exponential) and input parameters, appended
  experiment *E; // an experiment (data, initial conditions, inputs)
  experiment ref_E;
  int normalisation_type;
  gsl_matrix_int *norm_f; // (1|C) × F matrix
  gsl_matrix_int *norm_t; // (1|C) × F matrix
  gsl_vector *tmpF;
  gsl_vector *prior_tmp_a; // some memory allocation ..
  gsl_vector *prior_tmp_b; // .. for calculation buffers
  gsl_vector *prior_mu; // prior parameter: the medians of log normal distributions;
  gsl_matrix *prior_inverse_cov; // prior parameter: the inverse covariance in log space (where the prior is a Gaussian distribution);
  ode_solver *solver; // contains: cvode_mem; *odeModel; N_Vector y; N_Vector *yS; params;
} ode_model_parameters;

int ode_model_parameters_alloc(ode_model_parameters *omp, const problem_size *ps);
int ode_model_experiment_link(experiment *E);
int ode_model_parameters_free(ode_model_parameters *omp);
