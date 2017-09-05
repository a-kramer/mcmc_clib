#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "../mcmc/model_parameters_smmala.h"
#include "diagnosis_output.h"

int gsl_printf(const char *name, void *gsl_thing, int is_matrix){
  gsl_matrix *A;
  gsl_vector *x;
  int i,j,r,c;
  fflush(stdout);
  printf("[%s]",name);
  if (gsl_thing!=NULL){
    if (is_matrix){
      A=(gsl_matrix*) gsl_thing;
      r=(int) A->size1;
      c=(int) A->size2;
      printf(" %i × %i elements\n",r,c);
      for (i=0;i<r;i++){
	for (j=0;j<c;j++) printf("%g\t",gsl_matrix_get(A,i,j));
	printf("\n");
      }    
    }else{
      x=(gsl_vector*) gsl_thing;
      c=(int) x->size;
      printf(" %i elements\n",c);
      for (j=0;j<c;j++) printf("%g\t",gsl_vector_get(x,j));
      printf("\n");
    }
    printf("[/%s]\n",name);
  } else {
    printf("not set\n");
  }
  fflush(stdout);
  return GSL_SUCCESS;
}


int printf_omp(void *mp){
  int i,j,k,c,T,D,C;
  gsl_vector *t;
  char name[512];
  ode_model_parameters *omp;
  C=omp->size->C;
  D=omp->size->D;
  gsl_printf("p",omp->p,0);
  printf("Reference Experiment\n");
  gsl_printf("ref y init",omp->ref_E->init_y,0);
  t=omp->ref_E->t;
  if (t!=NULL){
    T=t->size;
    for (j=0; j<T; j++){
      printf("t[%i]\n",j);
      sprintf(name,"ref fy(t[%i])",j);
      //gsl_printf(sprintf("ref y(%i)",j),omp->ref_E->y[j],0);
      gsl_printf(name,omp->ref_E->fy[j],0);
      //gsl_printf(sprintf("ref yS(%i)",j),omp->ref_E->yS[j],1);
      sprintf(name,"ref fyS(t[%i])",j);
      gsl_printf(name,omp->ref_E->fyS[j],1);
    }
  }
  for (c=0;c<C;c++){
    printf("Experiment %i\n",c);
    gsl_printf("y_init",omp->E[c]->init_y,0);
    t=omp->E[c]->t;
    if (t!=NULL){
    T=t->size;
    for (j=0; j<T; j++){
      printf("t[%i]\n",j);
      //gsl_printf("y",omp->E[c]->y[j],0);
      sprintf(name,"E[%i] fy(t[%i])",c,j);
      gsl_printf(name,omp->E[c]->fy[j],0);
      sprintf(name,"E[%i] data(t[%i])",c,j);
      gsl_printf(name,omp->E[c]->data[j],0);
      sprintf(name,"E[%i] sd data(t[%i])",c,j);
      gsl_printf(name,omp->E[c]->sd_data[j],0);
      //gsl_printf("yS",omp->E[c]->yS[j],1);
      sprintf(name,"E[%i] fyS(t[%i])",c,j);
      gsl_printf(name,omp->E[c]->fyS[j],1);
      sprintf(name,"E[%i] oS(t[%i])",c,j);
      gsl_printf(name,omp->E[c]->oS[j],1);
    }
    } else {
      printf("E[%i] not properly defined; t missing\n",c);
    }    
  }
  fflush(stdout);
  return EXIT_SUCCESS;
}

