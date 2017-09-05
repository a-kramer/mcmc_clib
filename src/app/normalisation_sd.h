#include <gsl/gsl_matrix.h>
typedef struct {
  gsl_matrix *M;
  gsl_matrix *sd;
} gsl_matrix_sd;


int normalise_by_timepoint_with_sd(void *omp);
int normalise_by_state_var_with_sd(void *omp);
int ratio_with_sd(gsl_matrix_sd *A, gsl_matrix_sd *B);
