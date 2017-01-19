#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <regex.h>
#include <ctype.h>
#include <math.h>
#include "model_parameters_smmala.h"

// normalisation methods
#define DATA_IS_ABSOLUTE 0
#define DATA_NORMALISED_BY_REFERENCE 1
#define DATA_NORMALISED_BY_TIMEPOINT 2
#define DATA_NORMALISED_BY_STATE_VAR 3

// define target block types
#define INTEGER_BLOCK 1
#define DOUBLE_BLOCK 2

// data field ids, should be consecutive atm., because they are sometimes looped over.
#define i_time 0
#define i_reference_input 1
#define i_reference_data 2
#define i_sd_reference_data 3
#define i_input 4
#define i_data 5
#define i_sd_data 6
#define i_prior_mu 7
#define i_prior_icov 8
#define i_initial_conditions 9
#define i_norm_f 10
#define i_norm_t 11


typedef struct cnf_block_value_t cnf_block_value;
typedef struct cnf_block_value_t {
  value_t value;
  cnf_block_value *next;
};

typedef struct {
  char *library_file;
  char *output_file;
  double target_acceptance;
  double initial_stepsize;
  long sample_size;
  double t0;
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
  int id;
};



field_expression* field_expression_stack(int id,
					 field_expression *top,
					 regex_t *open,
					 regex_t *close);

field_expression* field_expression_init(field_names *fn);

int parse_config(FILE *cnf,
		 ode_model_parameters *omp,
		 problem_size *ps,
		 main_options *cnf_options);


int read_problem_definition(FILE *cnf,
			    ode_model_parameters *omp,
			    gsl_matrix_sd *RD,
			    const field_expression *fe,
			    const regex_t *comment,
			    problem_size *ps,
			    main_options *cnf_options);
