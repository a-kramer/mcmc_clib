#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <stdlib.h>
#include <stdio.h>
#include "../ode/ode_model.h"

// normalisation methods
#define DATA_IS_ABSOLUTE 0
#define DATA_NORMALISED_BY_REFERENCE 1
#define DATA_NORMALISED_BY_TIMEPOINT 2
#define DATA_NORMALISED_BY_STATE_VAR 3
#define DATA_NORMALISED_INDIVIDUALLY 4

// prior types
#define PRIOR_IS_UNKNOWN 0
#define PRIOR_IS_GAUSSIAN 1
#define PRIOR_IS_MULTIVARIATE 2
#define PRIOR_IS_GENERALISED 4
#define PRIOR_IS_KERNEL_DENSITY_ESTIMATE 8
#define PRIOR_IS_DVINE_COPULA 16
#define PRIOR_SIGMA_GIVEN 32
#define PRIOR_PRECISION_GIVEN 64
#define PTYPE(OPT,PROPERTY) (((OPT) & (PROPERTY)) > 0)

#define get_number_of_state_variables(omp) (omp->problem_size->N)
#define get_number_of_MCMC_variables(omp) (omp->problem_size->D)
#define get_number_of_model_parameters(omp) (omp->problem_size->P)
#define get_number_of_model_outputs(omp) (omp->problem_size->F)
#define get_number_of_model_inputs(omp) (omp->problem_size->U)
#define get_number_of_experimental_conditions(omp) (omp->problem_size->C)

typedef struct {
  int D; // number of sampling variables; unknown ODE-Model parameters
  int P; // number of total parameters of the ODE-Model (including dependent and known parameters)
  int N; // number of state varibales
  int F; // number of output functions
  int T; // number of measurement time-points for NORMALISATION_BY_REFERENCE case
  int U; // number of input parameters
  int C; // number of experimental conditions, excluding control (reference measurement)
} problem_size;

typedef struct {
  gsl_matrix *jacobian_y;
  gsl_matrix *jacobian_p;
  gsl_matrix *eJt;
  gsl_matrix *Jt;
  gsl_vector *tau;
  gsl_vector *x;
  gsl_vector *r;
} sensitivity_approximation; 

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
  gsl_matrix_view data_block_view; // view of this experiments data block inside the Data array, if that has been used
  gsl_matrix *data_block;  // either a pointer to a submatrix or self-allocated storage
  gsl_matrix_view sd_data_block_view; // standard deviation of this data block
  gsl_matrix *sd_data_block;
  gsl_vector_view *data_row; // for convenience, we define views for each row
  gsl_vector **data;        // data at time[j] is accessed as experiment[i]->data[j]
  gsl_vector_view *sd_data_row; // same for the standard deviations.
  gsl_vector **sd_data;
  gsl_vector **y; // y[t](i) vector of size T;
  gsl_vector **fy; // measurement model output functions
  gsl_vector **nfy; // normalising fy; temporary storage
  gsl_matrix **nfyS; // normalising fyS; temporary storage
  gsl_vector *init_y;
  gsl_vector *input_u;
  gsl_matrix *yS0;
  gsl_matrix **yS;
  gsl_matrix **fyS;
  gsl_matrix **oS; // observational sensitivities
  int NormaliseByExperiment;
  int NormaliseByTimePoint;
  gsl_vector_int *NormaliseByOutput;
} experiment;

typedef struct {
  gsl_vector *mu;
  gsl_vector **tmp; // temporary storage for prior calculations
  int n; // number of temporary vectors above
  union {
    gsl_matrix *Sigma_LU; // LU factors of the Sigma matrix
    gsl_matrix *inv_cov; // Same matrix already given in inverted form
    gsl_vector *sigma; // If Sigma is diagonal, this is the diagonal vector.
  };
  gsl_permutation *p; // used whenever LU decomposition is used
  int signum;         // also belongs to LU decomposition
  int type;           // for later use if maybe generalised prior is available, but also, diagonal/non-diagonal stuff;
} prior_t;

typedef struct {
  double t0;
  problem_size *size;
  gsl_vector *p; // memory for the ODE's parameters (they are
		 // exponential) and input parameters, appended: p=[exp(x),u];
  experiment **E; // an experiment (data, initial conditions, inputs, simulation results)
  experiment *ref_E;
  int normalisation_type;
  gsl_matrix_int *norm_f; // (1|C) × F matrix that stores normalisation information [old]
  gsl_matrix_int *norm_t; // (1|C) × F matrix that stores normalisation information [old]
  gsl_matrix *Data; // this is a block of size C*T*F, contains all data, if data is read as one huge array; otherwise this is NULL and data is stored in the experiment structure.
  gsl_matrix *sdData; // standard deviation of the datapoints
  gsl_vector *tmpF; // temporay storage of size F
  gsl_matrix *tmpDF; // temporary storage of size D×F
  prior_t *prior;
  gsl_matrix *FI_l; // likelihood contribution to fisher information matrix
  gsl_matrix *FI_p; // prior contribution to fisher information matrix
  sensitivity_approximation *S_approx; 
  ode_solver *solver; // contains: cvode_mem; *odeModel; N_Vector y; N_Vector *yS; params;
} ode_model_parameters;

int ode_model_parameters_alloc(ode_model_parameters *omp);
int ode_model_parameters_link(ode_model_parameters *omp);
int ode_model_parameters_free(ode_model_parameters *omp);
int init_E(ode_model_parameters *omp);
