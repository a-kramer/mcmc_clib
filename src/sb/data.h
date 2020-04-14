#ifndef DATA_H
#define DATA_H
#include <stdlib.h>
#include <glib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "norm.h"
#include <assert.h>

typedef struct data {
  int MajorIndex;
  int MinorIndex;
  int lflag;
  gsl_vector *input;
  gsl_vector *time;
  gsl_matrix *measurement;
  gsl_matrix *noise;
  gsl_vector *InitialValue;
  GString *Name;
} data_t;


void write_data(gpointer data, gpointer user_data);
GPtrArray* unwrap_data
(const GPtrArray *DataTable,
 const int *lflag,
 const experiment_type *ExperimentType,
 const GPtrArray *ExperimentName,
 const GPtrArray *input_override,
 const gsl_vector *default_time,
 const sbtab_t *Input,
 const sbtab_t *Output,
 map_t *IdxMap);

#endif
