#include <stdlib.h>
#include "h5block.h"

h5block_t* h5block_alloc(hsize_t rank){
  h5block_t *h5b=malloc(sizeof(h5block_t));
  h5b->size=calloc(rank,sizeof(hsize_t));
  h5b->chunk_size=calloc(rank,sizeof(hsize_t));
  h5b->offset=calloc(rank,sizeof(hsize_t));
  h5b->stride=calloc(rank,sizeof(hsize_t));
  h5b->count=calloc(rank,sizeof(hsize_t));
  h5b->block=calloc(rank,sizeof(hsize_t));
  return h5b;
}
void h5block_free(h5block_t *h5b){
  free(h5b->size);//calloc(rank,sizeof(hsize_t));
  free(h5b->chunk_size);//calloc(rank,sizeof(hsize_t));
  free(h5b->offset);//calloc(rank,sizeof(hsize_t));
  free(h5b->stride);//calloc(rank,sizeof(hsize_t));
  free(h5b->count);//calloc(rank,sizeof(hsize_t));
  free(h5b->block);//calloc(rank,sizeof(hsize_t));
}
