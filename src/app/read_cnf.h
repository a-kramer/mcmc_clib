#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <regex.h>
#ifdef _SMMALA
#include "model_parameters_smmala.h"
#endif
#ifdef _RMHMC
#include "model_parameters_rmhmc.h"
#endif

typedef struct {
  char *library_file;
  char *output_file;
  double target_acceptance;
  double initial_stepsize;
  long sample_size;
} main_options;

int ratio_with_sd(gsl_matrix *A, gsl_matrix *sdA, gsl_matrix *B, gsl_matrix *sdB);
int count_rows(FILE *cnf, regex_t *end, regex_t *comment);
int count_columns(const char *c);
int read_columns(char *c, double *vector, const int length);
int parse_config(FILE *cnf, ode_model_parameters *omp, main_options *cnf_options);
