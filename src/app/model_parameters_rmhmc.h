#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdlib.h>
#include <stdio.h>
#include "../ode/ode_model.h"

typedef struct {
  int D;
  gsl_vector *exp_x_u;  // memory for the exponential parameters and input parameters 
  double t0;
  gsl_vector *t;        // list of timepoints where measurements occur;
  gsl_matrix *Data;     //
  gsl_matrix *sdData;   // standard deviation of the datapoints
  int data_is_relative; // boolean: measurements are relative to a
			// reference measurement

  gsl_matrix *output_C; // linear output matrix y(t)=C*x(t); \dot x=f(t,x;\theta,u)
  gsl_vector *dl; // log likelihood gradient
  /* the variable cvodes_states contains the state vector returned by
   * the cvode integrator and needs memory to store y, yS, dyS, fy and fyS and dfyS;
   */
  double *cvodes_state; // i.e. N_Vector
  /* the variable cvodes_reference_states stores the same information
   * but for all time points 
   */
  double *cvodes_reference_state;
  gsl_vector *fy;
  gsl_matrix *fyS;
  gsl_matrix *dfyS;

  gsl_vector *rfy;
  gsl_matrix *rfyS;
  gsl_matrix *drfyS;

  gsl_matrix *oS;
  gsl_matrix *doS;

  gsl_matrix *input_u;  /* inputs applied in the lab: perturbations to
		       	 *  the model (double precision condition
			 *  parameters)
			 */
  gsl_vector *reference_u;
  
  gsl_vector *prior_tmp_a;       // some memory allocation ..
  gsl_vector *prior_tmp_b;       // .. for calculation buffers
  gsl_vector *prior_mu;          // prior parameter: the medians of log normal distributions;
  gsl_matrix *prior_inverse_cov; // prior parameter: the inverse covariance in log space;

  ode_solver *solver;        // contains: cvode_mem; *odeModel; N_Vector y; N_Vector *yS; params;
} ode_model_parameters;

int ode_model_parameters_alloc(ode_model_parameters *omp, int D, int N, int F,  int T, int U, int C);
int ode_model_parameters_free(ode_model_parameters *omp);