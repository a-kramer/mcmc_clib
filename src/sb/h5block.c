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


char* h5_to_char(hid_t g_id, const char *obj_name, const char *attr_name){
  int rank;
  hsize_t size[2]={1};
  herr_t status=0;
  char *str;
  H5T_class_t type_class;
  size_t type_size;
  /*attributes always have rank==1, at least when the H5LT API is used
    to make them*/
  h5_check(g_id,obj_name,attr_name);
  h5_get_info(g_id, obj_name, attr_name, &rank, size, &type_class, &type_size);
  if (attr_name){
    str = calloc(type_size,sizeof(char));
    status|=H5LTget_attribute_string(g_id,obj_name,attr_name,str);
  } else {
    str=calloc(type_size,sizeof(char));
    H5LTread_dataset_string(g_id, obj_name, str);
  }
  return str;
}


int h5_check(hid_t g_id, const char *obj_name, const char *attr_name){
  assert(obj_name);
  printf("[%s] checking existence of «%s».\n",__func__,obj_name);
  hid_t d_id=H5Dopen2(g_id, obj_name, H5P_DEFAULT);
  assert(H5LTfind_dataset(g_id, obj_name));
  if (attr_name){
    printf("[%s] checking existence of «%s».\n",__func__,attr_name);
    assert(H5LTfind_attribute(d_id, attr_name));
  }
  fflush(stdout);
  return 1;
}

void h5_get_info
(hid_t g_id,
 const char *obj_name,
 const char* attr_name,
 int *rank, hsize_t *size,
 H5T_class_t *type_class,
 size_t *type_size)
{
  herr_t status;
  assert(H5LTfind_dataset(g_id, obj_name));
  if (attr_name){
    rank[0]=1;
    status=H5LTget_attribute_info(g_id,obj_name,attr_name,size,type_class,type_size);
    assert(status>=0);
    fprintf(stderr,"[%s] attribute «%s» class %i, size %li rank %i (dims: %lli %lli).\n",
	   __func__,attr_name,*type_class, *type_size,
	   *rank,size[0],size[1]);
  } else {
    H5LTget_dataset_ndims(g_id, obj_name, rank);
    H5LTget_dataset_info(g_id, obj_name, size, type_class, type_size);
    fprintf(stderr,"[%s] dataset «%s»: class %i, size %li rank %i (dims: %lli %lli).\n",
	   __func__,obj_name,*type_class,*type_size,
	   *rank,size[0],size[1]);

  }
  
}

void* h5_to_gsl_int(hid_t g_id, const char *obj_name, const char* attr_name){
  hsize_t size[2]={1};
  herr_t status=0;
  gsl_vector_int *v;
  gsl_matrix_int *m;
  H5T_class_t type_class;
  size_t type_size;
  int rank;
  /*attributes always have rank==1, at least when the H5LT API is used
    to make them*/
  h5_check(g_id,obj_name,attr_name);
  h5_get_info(g_id, obj_name, attr_name, &rank, size, &type_class, &type_size);
  if (attr_name){
    v=gsl_vector_int_alloc(size[0]);
    status|=H5LTget_attribute_int(g_id,obj_name,attr_name,v->data);
    return v;
  } else {
    switch(rank){
    case 1:
      v=gsl_vector_int_alloc(size[0]);
      H5LTread_dataset_int(g_id, obj_name, v->data);
      return v;
    case 2:
      m=gsl_matrix_int_alloc(size[0],size[1]);
      H5LTread_dataset_int(g_id, obj_name, m->data);
      return m;
    default:
      fprintf(stderr,"[%s] currently, the number of dimensions should be either 1 or 2, got: %i.\n",__func__,rank);
      abort();      
    }
  }
}


void* h5_to_gsl(hid_t g_id, const char *obj_name, const char* attr_name){
  int rank;
  hsize_t size[2]={1};
  herr_t status=0;
  gsl_vector *v;
  gsl_matrix *m;
  H5T_class_t type_class;
  size_t type_size;
  /*attributes always have rank==1, at least when the H5LT API is used
    to make them*/
  h5_check(g_id,obj_name,attr_name);
  h5_get_info(g_id, obj_name, attr_name, &rank, size, &type_class, &type_size);
  if (attr_name){
    v=gsl_vector_alloc(size[0]);
    status|=H5LTget_attribute_double(g_id,obj_name,attr_name,v->data);
    return v;
  } else {
    switch (rank){
    case 1:
      v=gsl_vector_alloc(size[0]);
      H5LTread_dataset_double(g_id, obj_name, v->data);
      return v;
    case 2:
      m=gsl_matrix_alloc(size[0],size[1]);
      H5LTread_dataset_double(g_id, obj_name, m->data);
      return m;
    default:
      fprintf(stderr,"[%s] currently, the number of dimensions should be either 1 or 2, got: %i.\n",__func__,rank);
      abort();
    }
  }
  return v;
}

void gsl_matrix_to_h5(gsl_matrix *m, hid_t loc_id, const char *obj_name){
  herr_t status=0;
  assert(m);
  hsize_t size[2]={m->size1,m->size2};
  status=H5LTmake_dataset_double(loc_id, obj_name, 2, size, m->data);
  assert(status>0);
}


void gsl_vector_to_h5(gsl_vector *v, hid_t loc_id, const char *obj_name, const char *attr_name){
  herr_t status=0;
  hsize_t size=v->size;
  if (attr_name){
    status=H5LTset_attribute_double(loc_id, obj_name, attr_name, v->data, size);
  } else {
    status=H5LTmake_dataset_double(loc_id, obj_name, 1, &size, v->data);
  }
  assert(status>0);
}
