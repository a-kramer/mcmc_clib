#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <stdlib.h>
#include <stdio.h>
#include "../ode/ode_model.h"
// prior types
#include "ptype.h"
// normalisation methods
#define DATA_IS_ABSOLUTE 0
#define DATA_NORMALISED_BY_REFERENCE 1
#define DATA_NORMALISED_BY_TIMEPOINT 2
#define DATA_NORMALISED_BY_STATE_VAR 3
#define DATA_NORMALISED_INDIVIDUALLY 4
// for individual normalisation:
#define NEEDS_NORMALISATION(E) (((E)->NormaliseByExperiment>=0) || ((E)->NormaliseByTimePoint)>=0 || ((E)->NormaliseByOutput)!=NULL)

#define get_number_of_state_variables(omp) (omp->size->N)
#define get_number_of_MCMC_variables(omp) (omp->size->D)
#define get_number_of_model_parameters(omp) (omp->size->P)
#define get_number_of_model_outputs(omp) (omp->size->F)
#define get_number_of_model_inputs(omp) (omp->size->U)
#define get_number_of_experimental_conditions(omp) (omp->size->C)

#define set_number_of_state_variables(omp,N) (omp->size->N=N)
#define set_number_of_MCMC_variables(omp,D) (omp->size->D=D)
#define set_number_of_model_parameters(omp,P) (omp->size->P=P)
#define set_number_of_model_outputs(omp,F) (omp->size->F=F)
#define set_number_of_model_inputs(omp,U) (omp->size->U=U)
#define set_number_of_experimental_conditions(omp,C) (omp->size->C=C)


typedef struct {
  int D; // number of sampling variables; unknown ODE-Model parameters
  int P; // number of total parameters of the ODE-Model (including dependent and known parameters)
  int N; // number of state varibales
  int F; // number of output functions
  int T; // maximum size of time vectors over all experiments: T = max_c E[c]->t->size
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

/* helpful views of the data matrices, which can be used to address
 * rows as vectors. Once created, the views function in the background
 * and are not explicitely used.
 */
typedef struct {
  gsl_matrix_view data_block; // view of this experiments data block inside the Data array, if that has been used
  gsl_matrix_view sd_data_block; // standard deviation of this data block
  gsl_vector_view *data_row; // for convenience, we define views for each row
  gsl_vector_view *sd_data_row; // same for the standard deviations.
} view_t;


/* This structure contains some pre-allocated space for normalisation
 * purposes.
 */
typedef struct {
  gsl_vector **fy; // normalising fy; simulation result
  gsl_matrix **fyS; // normalising fyS; simulation result
  gsl_vector **data; // normalising data vector  (from a lab)
  gsl_vector **stdv; // and its standard deviation;  
} normalisation_t;

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
  view_t *view;
  normalisation_t *normalise;
  int lflag;                // TRUE if this data should influence LogLikelihood explicitely; FALSE if it is used in some implicit way (e.g. normalisation) 
  gsl_matrix *data_block;   // either a pointer to a submatrix or self-allocated storage
  gsl_matrix *sd_data_block;// either a pointer to a submatrix or self-allocated storage
  gsl_vector **data;        // data at time[j] is accessed as experiment[i]->data[j]
  gsl_vector **sd_data;
  gsl_vector *p; // working memory, to store this experiments parameter vector
  gsl_vector **y; // y[t](i) vector of size T;
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
} experiment;

typedef struct {
  union{
    gsl_vector *mu;
    gsl_vector *alpha;
  };
  gsl_vector **tmp; // temporary storage for prior calculations
  int n; // number of temporary vectors above
  union {
    gsl_matrix *Sigma_LU; // LU factors of the Sigma matrix
    gsl_matrix *inv_cov; // Same matrix already given in inverted form
    gsl_vector *sigma; // If Sigma is diagonal, this is the diagonal vector.
    gsl_vector *beta;  // this is a parameter of generalised Gaussians
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
  int normalisation_type;
  gsl_matrix *Data; // this is a block of size C*T*F, contains all data, if data is read as one huge array; otherwise this is NULL and data is stored in the experiment structure directly.
  gsl_matrix *sdData; // standard deviation of the datapoints above, NULL if Data is NULL;
  gsl_vector *tmpF; // temporay storage of size F
  gsl_matrix *tmpDF; // temporary storage of size D×F
  prior_t *prior;
  sensitivity_approximation **S_approx; 
  ode_solver **solver; // contains: cvode_mem; *odeModel; N_Vector y; N_Vector *yS; params; per experiment
} ode_model_parameters;

int ode_model_parameters_alloc(ode_model_parameters *omp);
int ode_model_parameters_link(ode_model_parameters *omp);
int ode_model_parameters_free(ode_model_parameters *omp);
int init_E(ode_model_parameters *omp);
