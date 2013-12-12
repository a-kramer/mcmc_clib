#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <regex.h>
#include "model_parameters2.h"

int ratio_with_sd(gsl_matrix *A, gsl_matrix *sdA, gsl_matrix *B, gsl_matrix *sdB);
int count_columns(const char *c);
int read_columns(char *c, double *vector, const int length);
int parse_config(FILE *cnf, ode_model_parameters *omp);
