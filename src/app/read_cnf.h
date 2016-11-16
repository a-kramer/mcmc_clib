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
} main_options;  // these are user supplied options to the program

typedef struct  {
  int n; //number of field names
  size_t max_length; // maximum length of fields
  char **name;
} field_names;
int field_names_init(field_names *fn);

typedef struct field_expression_t field_expression;
struct field_expression_t {
  regex_t *opening_bracket; // regular expression for [field_name]
  regex_t *closing_bracket; // regular expression for [/field_name]
  field_expression *next;
};

field_expression* field_expression_stack(field_expression *top, regex_t *open, regex_t *close);
field_expression* field_expression_init(field_names *fn);

int ratio_with_sd(gsl_matrix *A, gsl_matrix *sdA, gsl_matrix *B, gsl_matrix *sdB);
int count_rc(FILE *cnf, regex_t *end, regex_t *comment, int *rows, int *columns);
int count_columns(const char *c);
int read_vector(char *c, double *vector, const int length);
int parse_config(FILE *cnf, ode_model_parameters *omp, main_options *cnf_options);
int determine_problem_size(FILE *cnf, const field_expressions *fe, const regex_t *comment, problem_size *ps, main_options *cnf_options);
int read_problem_definition(FILE *cnf, const field_expressions *fe, const regex_t *comment, problem_size *ps, main_options *cnf_options);
