#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mpi.h>
#include "../mcmc/model_parameters_smmala.h"
#include "diagnosis_output.h"
#include "../mcmc/ptype.h"
// type is an integer of bitwise flags
// it can contain multiple simultaneously activated flags
// so, type can be true for property1 and ALSO property2 among other things

typedef union {
  double d;
  int i;
} val;

typedef union {
  gsl_vector *d;
  gsl_vector_int *i;  
} gsl_vec;

typedef union {
  gsl_matrix *d;
  gsl_matrix_int *i;  
} gsl_mat;


int gsl_fprintf(FILE *a,const char *name, void *gsl_thing, int type){
  gsl_mat A;
  gsl_vec x;
  int i,j,r,c;
  val v;
  //assert(name);
  fprintf(a,"[%s]",name);
  if (gsl_thing){
    if (PTYPE(type,GSL_IS_MATRIX)){
      if (PTYPE(type,GSL_IS_DOUBLE)){
	A.d=(gsl_matrix*) gsl_thing;
	r=(int) A.d->size1;
	c=(int) A.d->size2;
      } else if (PTYPE(type,GSL_IS_INT)){
	A.i=(gsl_matrix_int*) gsl_thing;
	r=(int) A.i->size1;
	c=(int) A.i->size2;
      } else {
	fprintf(stderr,"[gsl_fprintf] GSL matrix has unknown element type: «%x»",type);
	return GSL_EINVAL;
      }
      fprintf(a," %i × %i elements\n",r,c);
      if (PTYPE(type,GSL_IS_DOUBLE)){
	for (i=0;i<r;i++){
	  for (j=0;j<c;j++) {	    
	    v.d=gsl_matrix_get(A.d,i,j);
	    fprintf(a,"%+g\t",v.d);
	  }
	  fprintf(a,"\n");
	}
      } else if (PTYPE(type,GSL_IS_INT)){
	for (i=0;i<r;i++){
	  for (j=0;j<c;j++) {	    
	    v.i=gsl_matrix_int_get(A.i,i,j);
	    fprintf(a,"%i\t",v.i);
	  }
	  fprintf(a,"\n");
	}
      } else {
	fprintf(stderr,"[gsl_fprintf] GSL matrix has unknown element type: «%x»",type);
	return GSL_EINVAL;
      }
    } else if (PTYPE(type,GSL_IS_VECTOR)){
      if (PTYPE(type,GSL_IS_DOUBLE))	{
	x.d=(gsl_vector*) gsl_thing;
	c=(int) x.d->size;
      } else if (PTYPE(type,GSL_IS_INT)) {
	x.i=(gsl_vector_int*) gsl_thing;
	c=(int) x.i->size;
      }
      fprintf(a," %i elements\n",c);
      if (PTYPE(type,GSL_IS_DOUBLE)){
	for (j=0;j<c;j++){
	  v.d=gsl_vector_get(x.d,j);
	  fprintf(a,"%+g\t",v.d);
	}
	fprintf(a,"\n");
      }else if (PTYPE(type,GSL_IS_INT)){
	for (j=0;j<c;j++){
	  v.i=gsl_vector_int_get(x.i,j);
	  fprintf(a,"%i\t",v.i);
	}
	fprintf(a,"\n");
      } else {
	fprintf(stderr,"[gsl_fprintf] GSL vector has unknown element type: «%x»",type);
	return GSL_EINVAL;
      }
    } else {
      fprintf(stderr,"[gsl_fprintf] unknown GSL entity type: «%x»",type);
      return GSL_EINVAL;
    }
    fprintf(a,"[/%s]\n",name);
  } else {
    fprintf(a,"[gsl_fprintf] GSL pointer is not set\n");
    return GSL_EINVAL;
  }
  return GSL_SUCCESS;
}

int gsl_printf(const char *name, void *gsl_thing, int type){
  return gsl_fprintf(stdout,name,gsl_thing,type);
}

int printf_omp(void *mp){
  int j,c,T,C;
  gsl_vector *t;
  char name[512];
  ode_model_parameters *omp=mp;
  FILE *f;
  int r;
  MPI_Comm_rank(MPI_COMM_WORLD,&r);
  sprintf(name,"ode_model_state_of_rank_%i.txt",r);
  f=fopen(name,"w");
  C=omp->size->C;
  fprintf(f,"[printf_omp]\n"); fflush(stdout);
  gsl_fprintf(f,"p",omp->p,GSL_IS_DOUBLE | GSL_IS_VECTOR); fflush(stdout);
  for (c=0;c<C;c++){
    fprintf(f,"Experiment %i\n",c);
    gsl_fprintf(f,"y_init",omp->E[c]->init_y,GSL_IS_DOUBLE | GSL_IS_VECTOR);
    t=omp->E[c]->t;
    if (t!=NULL){
      T=t->size;
      for (j=0; j<T; j++){
	fprintf(f,"t[%i]\n",j);
	gsl_fprintf(f,"y",omp->E[c]->y[j],GSL_IS_DOUBLE | GSL_IS_VECTOR);
	sprintf(name,"E[%i] fy(t[%i])",c,j);
	gsl_fprintf(f,name,omp->E[c]->fy[j],GSL_IS_DOUBLE | GSL_IS_VECTOR);
	sprintf(name,"E[%i] data(t[%i])",c,j);
	gsl_fprintf(f,name,omp->E[c]->data[j],GSL_IS_DOUBLE | GSL_IS_VECTOR);
	sprintf(name,"E[%i] sd data(t[%i])",c,j);
	gsl_fprintf(f,name,omp->E[c]->sd_data[j],GSL_IS_DOUBLE | GSL_IS_VECTOR);
	sprintf(name,"E[%i] fyS(t[%i])",c,j);
	gsl_fprintf(f,name,omp->E[c]->fyS[j],GSL_IS_DOUBLE | GSL_IS_MATRIX);
	sprintf(name,"E[%i] oS(t[%i])",c,j);
	gsl_fprintf(f,name,omp->E[c]->oS[j],GSL_IS_DOUBLE | GSL_IS_MATRIX);
      }
    } else {
      fprintf(f,"E[%i] not properly defined; t missing\n",c);
    }    
  }
  fclose(f);
  return EXIT_SUCCESS;
}

