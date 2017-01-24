#include "dynamic_array.h"

var_ndim_array* nda_alloc(){
  var_ndim_array *new_nda;
  new_nda=(var_ndim_array *)malloc(sizeof(var_ndim_array));
  new_nda->root=NULL;
  new_nda->N=0;
  new_nda->nd=0;
  new_nda->size=NULL;
  return new_nda;
}

int nda_free(var_dim_array *nda){
  if (nda->root == NULL){
    free(nda->size);
    free(nda);
  } else {
    perror("var dim array is not empty.");
    exit(-1);
  }
  return GSL_SUCCESS;
}

int nda_push(var_ndim_array* nda, value_t* value){
  cnf_block_value *p;
  p=nda->root;
  
  p=(cnf_block_value*) malloc(sizeof(cnf_block_value));
  if (nda->btype==DOUBLE_BLOCK){
    p->value->d=value->d;
  }else if (nda->btype==INTEGER_BLOCK){
    p->value->i=value->i;
  } else {
    fprintf(stderr,"nda_push: error for target type: %i\n",nda->btype);
    exit(-3);
  }
  nda->N++;
  p->next=nda->root;
  nda->root=p;
  return GSL_SUCCESS;
}
int nda_pop(var_ndim_array* nda, value_t* value){
  cnf_block_value *p;
  p=nda->root;
  if (p!=NULL){
    if (nda->btype==DOUBLE_BLOCK){
      value->d=p->value->d;
    }else if (nda->btype==INTEGER_BLOCK){
      value->i=p->value->i;
    } else {
      fprintf(stderr,"nda_pop: error for target type: %i\n",nda->btype);
      exit(-3);
    }
  }else{
    return GSL_EINVAL;
  }
  nda->root=p->next;
  free(p);
  nda->N--;
  return GSL_SUCCESS;
}




/* nda_to_gsl gets a variable length n-dim array and copies the
 * contents, which it allocates itself.  This matrix directly
 * represents the contents of a file. The data can later be copied
 * into other structures. The return object will always be a matrix,
 * because some cfg blocks may be a matrix in nature (a table), but
 * because C=1 and T=1 have only one row but are expected to be a
 * matrix nontheless.
 */
int nda_to_gsl(var_dim_array *nda, gsl_object *G){
  int i,N=nda->N;
  int nd=nda->nd;
  int r,c;
  //  gsl_object *G;

  value_t value;
  if (nda->btype==DOUBLE_BLOCK){
    double *data;
    G->gsl->is_double++;
    r=nda->size[0];
    c=nda->size[1];
    G->gsl->matrix=gsl_matrix_alloc(r,c);
    G->gsl->is_matrix++;
    data=G->gsl->matrix->data;
    for (i=1;i<=N;i++) {
      if (nda_pop(nda,value)!= GSL_EINVAL){
	data[N-i]=value.d;
      } else {
	fprintf(stderr,"[nda_to_gsl] array contained fewer values than expected\n");
	exit(-3);	
      }
    }
  } else if (nda->btype==INTEGER_BLOCK){
    int *data;
    G->gsl->is_int++;
    r=nda->size[0];
    c=nda->size[1];
    G->gsl->matrix_int=gsl_matrix_int_alloc(r,c);
    G->is_matrix++;
    data=G->gsl->matrix_int->data;
    for (i=1;i<=N;i++) {
      if (nda_pop(nda,value)!= GSL_EINVAL){
	data[N-i]=value.i;
      } else {
	perror("[nda_to_gsl] array contained fewer values than expected");
	exit(-3);	
      }
    }
  } else {
    perror("nda_to_gsl: unknown block type: %i\n",nda->btype);
    exit(-3);    
  }
  return GSL_SUCCESS;
}

/* Same as nda_to_gsl but checks whether nda represents a matrix or a
 * vector. The gsl object will then either be a gsl_vector or a
 * gsl_matrix (possibly with _int suffix).  This can be used, whenever
 * it is important to make such a distinction: i.e. whether a block
 * contains one line or several makes a qualitative difference to the
 * setup.
 */

int nda_to_gsl_a(var_dim_array *nda, gsl_object *G){
  int i,N=nda->N;
  int nd=nda->nd;
  int r,c;
  //  gsl_object *G;

  value_t value;
  if (nda->btype==DOUBLE_BLOCK){
    double *data;
    G->gsl->is_double++;
    if (nd==2) {
      r=nda->size[0];
      c=nda->size[1];
      G->gsl->matrix=gsl_matrix_alloc(r,c);
      G->gsl->is_matrix++;
      data=G->gsl->matrix->data;
    } else if (nd==1) {
      r=nda->size[0];
      G->gsl->vector=gsl_vector_alloc(r);
      G->gsl->is_vector++;
      data=G->gsl->vector->data;
    } else {
      fprintf(stderr,"[nda_to_gsl] cannot tell if data is vector or matrix: nd=%i\n",nda->nd);
      exit(-3);    
    }
    for (i=1;i<=N;i++) {
      if (nda_pop(nda,value)!= GSL_EINVAL){
	data[N-i]=value.d;
      } else {
	fprintf(stderr,"[nda_to_gsl] array contained fewer values than expected\n");
	exit(-3);	
      }
    }
  } else if (nda->btype==INTEGER_BLOCK){
    int *data;
    G->gsl->is_int++;
    if (nd==2) {
      r=nda->size[0];
      c=nda->size[1];
      G->gsl->matrix_int=gsl_matrix_int_alloc(r,c);
      G->is_matrix++;
      data=G->gsl->matrix_int->data;
    } else if (nd==1) {
      r=nda->size[0];
      G->gsl->vector_int=gsl_vector_int_alloc(r);
      G->gsl->is_vector++;
      data=G->gsl->vector_int->data;
    } else {
      fprintf(stderr,"[nda_to_gsl] cannot tell if int_data is vector or matrix: nd=%i\n",nda->nd);
      exit(-3);    
    }
    for (i=1;i<=N;i++) {
      if (nda_pop(nda,value)!= GSL_EINVAL){
	data[N-i]=value.i;
      } else {
	perror("[nda_to_gsl] array contained fewer values than expected");
      exit(-3);	
      }
    }
  } else {
    perror("nda_to_gsl: unknown block type: %i\n",nda->btype);
    exit(-3);    
  }
  return GSL_SUCCESS;
}
