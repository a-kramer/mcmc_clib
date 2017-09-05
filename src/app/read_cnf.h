#include <regex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "../mcmc/model_parameters_smmala.h"
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
#define i_ref_initial_conditions 10
#define i_norm_f 11
#define i_norm_t 12
#define NumberOfFields 13

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

int read_problem_definition(FILE *cnf,
			    ode_model_parameters *omp,
			    const field_expression *fe,
			    const regex_t *comment,
			    main_options *cnf_options);

int parse_config(FILE *cnf,
		 ode_model_parameters *omp,
		 main_options *cnf_options);
