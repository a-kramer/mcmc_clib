#ifndef NORMALISATION_H
#define NORMALISATION_H
#include <stdlib.h>
#include <glib.h>
#include "sbtab.h"
#include "h5block.h"

#define UNMAPPED -1
#define ISMAPPED(i) ((i)>=0)

typedef enum e_type {unknown_type,dose_response,time_series} experiment_type;
typedef GArray* index_map_t;

typedef struct norm {
  GArray* output;
  GArray* time;
  GArray* experiment;
} norm_t;

typedef struct ExperimentIndexMap {
  GArray *major;
  GArray *minor;
  GPtrArray *flat;
} map_t;

typedef struct data_attribute_buffer {
  h5block_t *h5data;
  h5block_t *h5stdv;
  map_t *IdxMap;
  norm_t *N;
} data_attr_buffer_t;

map_t* empty_map(int nE);
void map_free(map_t *m);
int major_index(const map_t *m, int flat_index);
int minor_index(const map_t *m, int flat_index);
int flat_index(const map_t *m, int major, int minor);
void map_index(map_t *m, int major, int minor);


norm_t* norm_alloc(int nE, int nO);
int norm_index(index_map_t M, int i, int *out_j);
void norm_index_append(index_map_t M, int j);
norm_t* normalisation(sbtab_t *ExperimentTable, experiment_type *ExperimentType, sbtab_t *OutputTable, GPtrArray *DataTable, map_t *IdxMap);
#endif
