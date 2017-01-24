#include "read_cnf.h"

typedef struct {
  int N; //number of elements; always correct, when using push and pop
  int nd; // number of dimensions; 
  int *size; // size of each dimension;these are set when a data block is
	  // read from file
  int btype;
  cnf_block_value *root; // root element for the variable length array
} var_ndim_array;


typedef union{
  double d;
  int i;
} value_t;

// one gsl_object which can be a matrix or vector, double or integer
// this is used while reading data from a file and storing it for a bit.
typedef union {
  gsl_matrix *matrix;
  gsl_vector *vector;
  gsl_matrix_int *matrix_int;
  gsl_vector_int *vector_int;
} gsl_t;

typedef struct {
  int is_double;
  int is_int;
  int is_matrix;
  int is_vector;
  gsl_t *gsl;
} gsl_object;

var_ndim_array* nda_alloc();
int nda_push(var_ndim_array* nda, value_t* value);
int nda_to_gsl(var_dim_array *nda, gsl_object *G);
int nda_to_gsl_a(var_dim_array *nda, gsl_object *G);
