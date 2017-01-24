#include "model_parameters_smmala.h"

typedef struct {
  gsl_matrix *M;
  gsl_matrix *sd;
} gsl_matrix_sd;


int normalise_by_timepoint_with_sd(ode_model_parameters *omp, problem_size *ps);
int normalise_by_state_var_with_sd(ode_model_parameters *omp, problem_size *ps);
int ratio_with_sd(gsl_matrix_sd *A, gsl_matrix_sd *B);
