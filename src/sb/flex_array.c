#include "flex_array.h"
#include <stdio.h>

flex_double* flex_array_alloc(size_t n){
  flex_double *a;
  a=malloc(sizeof(flex_double));
  a->data=malloc(sizeof(double)*n);
  a->length=n;
  a->numel=0;
  return a;
}

int flex_array_resize(flex_double *a, size_t n){
  int status=EXIT_SUCCESS;
  if (a!=NULL && a->data!=NULL) {
    a->data=realloc(a->data,sizeof(double)*n);
    if (a->data==NULL) status=EXIT_FAILURE;
    a->length=n;
  }
  return status;
}

int flex_array_set(flex_double *a, size_t i, double d){
  int status=EXIT_SUCCESS;
  size_t m;
  if (a!=NULL){
    if (a->length<i){
      m=a->length + a->length/10;
      m=(m<i)?i+1:m;
      status&=flex_array_resize(a,m);
    }
    if (status!=EXIT_FAILURE) {
      a->data[i]=d;
    }
    if (i>=a->numel) a->numel++;
  } else {
    status=EXIT_FAILURE;
  }
  return status;
}

double flex_array_get(flex_double *a, size_t i){
  double d;
  if (a!=NULL && a->data!=NULL && i<a->numel){
    d=a->data[i];
  } else {
    fprintf(stderr,"[flex_array_get] ");
    if (a==NULL) fprintf(stderr,"array is not initialised.\n");
    else if (a->data==NULL) fprintf(stderr,"array is broken, a->data==NULL\n");
    else fprintf(stderr,"size(a)=[%i of %i]; getting a[%i]; size violation. exit.\n",a->numel,a->length,i);
    exit(-1);
  }
  return d;
}



int flex_array_free(flex_double *a){
  if (a!=NULL){
    if (a->data!=NULL) free(a->data);
    free(a->);
  }
  return EXIT_SUCCESS;
}
