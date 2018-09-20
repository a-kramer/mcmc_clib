#include <gsl/gsl_matrix.h>
typedef struct {
  gsl_matrix *M;
  gsl_matrix *sd;
} gsl_matrix_sd;


int normalise_with_sd(void *model_parameters);
