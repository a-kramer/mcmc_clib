#ifndef DATA_H
#define DATA_H
#include <glib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "sbtab.h"
// Experiment Types
typedef enum e_type {unknown_type,dose_response,time_series} experiment_type;
typedef struct {
  GArray *output;
  Garray *time;
  GArray *experiment;
} normalisation_t;

typedef struct data {
  int MajorIndex;
  int MinorIndex;
  int lflag;
  gsl_vector *input;
  gsl_vector *time;
  int NormaliseByExperiment;
  gsl_vector_int NormaliseByOutput;
  gsl_matrix *measurement;
  gsl_matrix *noise;
  gsl_vector *InitialValue;
} data_t;




#endif
