#ifndef H5BLOCK_H
#define H5BLOCK_H
#include "hdf5.h"
#include "hdf5_hl.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/* collects all ids and arrays the hdf5 needs
 */
typedef struct {
	hid_t file_id;
	hid_t group_id;
	hid_t dataset_id;
	hid_t chunk_id;
	hid_t dataspace_id;
	hsize_t *size;
	hsize_t *chunk_size;
	hsize_t *offset;
	hsize_t *stride;
	hsize_t *count;
	hsize_t *block;
} h5block_t;

h5block_t* h5block_alloc(hsize_t rank);
void h5block_free(h5block_t *h5b);

char* h5_to_char(hid_t g_id, const char *obj_name, const char *attr_name);
int h5_check(hid_t g_id, const char *obj_name, const char *attr_name);
void h5_get_info(hid_t g_id, const char *obj_name, const char* attr_name, int *rank, hsize_t *size, H5T_class_t *type_class, size_t *type_size);
void* h5_to_gsl_int(hid_t g_id, const char *obj_name, const char* attr_name);
void* h5_to_gsl(hid_t g_id, const char *obj_name, const char* attr_name);
void gsl_matrix_to_h5(gsl_matrix *m, hid_t loc_id, const char *obj_name);
void gsl_vector_to_h5(gsl_vector *v, hid_t loc_id, const char *obj_name, const char *attr_name);

#endif
