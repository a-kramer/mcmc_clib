#ifndef H5BLOCK_H
#define H5BLOCK_H
#include "hdf5.h"
#include "hdf5_hl.h"

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
#endif
