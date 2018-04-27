#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_int.h>
#include "normalisation_sd.h"
#include "../mcmc/model_parameters_smmala.h"
/* All normalisation functions in this file use absolte values of
 * derivatives and add them weighted by individual sds. This works
 * under the assumption that the result's standard deviation is more
 * accurately estimated by: sd(r=d/s) = |r/d|*sd(d) + |r/s|*sd(s).
 * This overestimates uncorrelated errors slightly compared to
 * (sqrt((r/d)²*sd(d))² + (r/s)²*sd(s)²)), which is the Taylor series based
 * expression.  However, in systems biology, the measurement errors
 * are often large, so the bias in a truncated Taylor series is
 * significant. Therefore, we err on the side of caution.  We also assume
 * both d and s to be positive, hence absolute values by default.
 */

int normalise_by_timepoint_with_sd(void *mp){
  /* here normalisation is just one entry, containing the
   * normalisation time index. Many of the operation are done in
   * place, so variable interpretation changes.
   */
  int r,c,j,l;
  ode_model_parameters *omp=mp;
  int C=omp->size->C;
  int T=omp->size->T;
  int F=omp->size->F;
  gsl_vector *d,*sd_d,*sd_s,*s;
  gsl_vector *tmp; // intermediate results

  tmp=omp->tmpF; //gsl_vector_alloc(F);
 
  r=omp->norm_t->size1;
  /*printf("r=%i\nnorm_t=[\n",r);
  for (c=0;c<r;c++){
    for (j=0;j<(omp->norm_t->size2);j++) printf("%i",gsl_matrix_int_get(omp->norm_t,c,j));
    printf("\n");
  }
  printf("]\n");
  */
  for (c=0;c<C;c++){
    l=gsl_matrix_int_get(omp->norm_t,c%r,0);
    T=omp->E[c]->t->size;
    //printf("l=%i\n",l);
    fflush(stdout);
    for (j=0;j<T;j++) {
      d=omp->E[c]->data[j];
      sd_d=omp->E[c]->sd_data[j];
      s=omp->E[c]->data[l];
      sd_s=omp->E[c]->sd_data[l];
      gsl_vector_div(d,s);
      // calculate standard deviation of result:
      gsl_vector_set_zero(tmp);
      gsl_vector_add(tmp,d);
      // add both error terms
      gsl_vector_div(sd_s,s);
      gsl_vector_mul(tmp,sd_s); // so: (d/s) * (sd_s/s) = (d/s²) sd_s   
        
      gsl_vector_div(sd_d,s);
      gsl_vector_add(sd_d,tmp); // so: (sd_d/s) + sd_s × (d/s²) = (d/s)×(sd_d/d) + (d/s)×(sd_s/s);
    }
  }
  /* all datapoints are now normalised in place. Some data points are
   * now exactly 1, but the simulations will be as well, so they won't
   * contribute to the likelihood.
   */
  //  gsl_vector_free(tmp);
  return GSL_SUCCESS;  
}

int normalise_by_state_var_with_sd(void *mp){
  /* here normalisation consists of two lines. Line 1 lists the state
   * variables to use for normalisation. Line 2 selects the time index
   * of that state variable to normalise at.
   */
  int c,i,j;
  int rows_f, rows_t;
  gsl_vector *tmp; // intermediate results
  //gsl_vector_view D,SD;
  gsl_vector *r,*sd_r,*d,*sd_d;
  ode_model_parameters *omp=mp;
  int i_f, i_t;
  int C=omp->size->C;
  int T=omp->size->T;
  int F=omp->size->F;
  // here, we need to allocate some memory for s and sd_s, because the
  // normalisation is more complex and different from state variable
  // to state variable.
  tmp=omp->tmpF;
  r=omp->ref_E->fy[0];
  sd_r=gsl_vector_alloc(F);
  rows_t=omp->norm_t->size1;
  rows_f=omp->norm_f->size1;
  printf("normalising with state and time index (rows_f=%i, rows_t=%i)\n",rows_f,rows_t); fflush(stdout);
  //printf("\n");
  //printf("# norm_f:\n");
  //gsl_matrix_int_fprintf(stdout,omp->norm_f,"%i");
  //printf("# norm_t:\n");
  //gsl_matrix_int_fprintf(stdout,omp->norm_t,"%i");
  
  for (c=0;c<C;c++) {
    T=omp->E[c]->t->size;
    //printf("# [norm_f/norm_t] c=%i/%i\n",c,C);
    // s is defined per block of experimental conditions (c=const.)
    for (i=0;i<F;i++) {
      //printf("# [norm_f/norm_t] i=%i/%i\n",i,F);

      // fill s with the appropriate data point:
      i_f=gsl_matrix_int_get(omp->norm_f,c%rows_f,i); 
      // at time index i_t, as specified by the normalisation matrix
      i_t=gsl_matrix_int_get(omp->norm_t,c%rows_t,i);
      //printf("i_t=%i\ti_f=%i\n",i_t,i_f);fflush(stdout);
      if (i_f<F && i_t<T){
	d=omp->E[c]->data[i_t];
	sd_d=omp->E[c]->sd_data[i_t];
	gsl_vector_set(r,i,gsl_vector_get(d,i_f));
	gsl_vector_set(sd_r,i,gsl_vector_get(sd_d,i_f));
      } else {
	printf("normalisation index pair out of bounds: 0≤%i<%i and 0≤%i<%i.\n",i_f,F,i_t,T); exit(-1);
      }
    }
    for (j=0;j<T;j++) { 
      d=omp->E[c]->data[j];
      sd_d=omp->E[c]->sd_data[j];
      // data will be scaled by 1/s
      gsl_vector_div(d,r);
      // calculate standard deviation of result:
      gsl_vector_memcpy(tmp,d);
      // add both error terms
      gsl_vector_div(sd_r,r);
      gsl_vector_mul(tmp,sd_r);
      // so: (d/s) * (sd_r/r) = (d/r²) sd_r   
      gsl_vector_div(sd_d,r);
      gsl_vector_add(sd_d,tmp);
      // so: (sd_d/r) + sd_r × (d/r²) = (d/r)×(sd_d/d) + (d/r)×(sd_r/r);
    }
  }
  gsl_vector_free(sd_r);
  return GSL_SUCCESS;  
}

/* we take the ratio A./B (by elements) B is repeated to match the
 * size of A. Both A and B contain a struct member sd, its standard
 * deviation values. A must have the same number of columns as B and
 * an integer multiple of rows, e.g. A=[B;B;B] (GNU Octave
 * Notation).  the result is stored in A.  All matrices are expected
 * to contain only positive elements.
 */
int ratio_with_sd(gsl_matrix_sd *A, gsl_matrix_sd *B){
  int l;
  int C,T,F;
  gsl_matrix_view a;
  gsl_matrix *R; // for the relative error of B

  T=B->M->size1;
  C=(A->M->size1)/T;
  F=A->M->size2;
  //  printf("sizes: T=%i, C=%i, F=%i\n",T,C,F);
  R=gsl_matrix_alloc(C*T,F);
  
  
  for (l=0;l<C;l++){
    a=gsl_matrix_submatrix(A->M,l*T,0,T,F);
    gsl_matrix_div_elements(&(a.matrix),B->M); //A'=A./B
  }
  
  for (l=0;l<C;l++){
    a=gsl_matrix_submatrix(A->sd,l*T,0,T,F);//sd(A)'=sd(A)./B
    gsl_matrix_div_elements(&(a.matrix),B->M);
  }
  
  gsl_matrix_div_elements(B->sd,B->M); //sd(B)'=sd(B)./B
  gsl_matrix_memcpy(R,A->M); // copy: R=A'=A./B
  for (l=0;l<C;l++){
    a=gsl_matrix_submatrix(R,l*T,0,T,F);
    gsl_matrix_mul_elements(&(a.matrix),B->sd); //R ends up being A./B*sd(B)/B
  }
  gsl_matrix_add(A->sd,R); // sd(A)./B + A./B*sd(B)/B
  gsl_matrix_free(R);
  return EXIT_SUCCESS;
}

int get_normalising_vector_with_sd(experiment *E, experiment *ref_E, gsl_vector **ref_data, gsl_vector **ref_stdv){
  //int i=E->NormaliseByExperiment;
  int t;
  int f,j,nj,k,nk;
  double val, dval;
  gsl_vector *v, *dv;
  assert(ref_data!=NULL && ref_stdv!=NULL);
  j=E->NormaliseByTimePoint;
  
  if(ref_E==NULL){
    ref_E=E;
    assert(j>=0);
  }
  int F=E->data[0]->size;
  int T=E->t->size;
  int rT=ref_E->t->size;

  assert(j<T);
  
  // v and dv are used to handle the NormaliseByTimePoint property
  for (t=0;t<T;t++){
    // find the right reference data and stdv
    if (j>=0){
      v=ref_E->data[j];
      dv=ref_E->sd_data[j];
    } else {
      v=ref_E->data[t%T];
      dv=ref_E->sd_data[t%T];
    }  
    if (E->NormaliseByOutput==NULL){
      // just link fy and fyS with the reference objects:
      ref_data[t]=v;
      ref_stdv[t]=dv;
    }else{ // each output is normalised differently
      ref_data[t]=E->normalise->data[t];
      ref_stdv[t]=E->normalise->stdv[t];
      // copy elements to fy and fyS
      for (f=0;f<F;f++){
	// output function:
	k=gsl_vector_int_get(E->NormaliseByOutput,f);
	val=gsl_vector_get(v,k);
	dval=gsl_vector_get(dv,k);
	gsl_vector_set(ref_data[t],f,val);
	gsl_vector_set(ref_stdv[t],f,val);
      }
    }
  }
  return GSL_SUCCESS;
}

int normalise_with_sd(void *model_parameters){
  ode_model_parameters *mp=model_parameters;
  int c,j,k,l;
  int C=mp->size->C;
  int T=mp->size->T;
  int D=mp->size->D;
  int F=mp->size->F;
  assert(T>0);
  gsl_vector **data, **stdv;
  gsl_vector *r_data[T], *r_stdv[T];
  gsl_vector *v;
  int i,nt;
  experiment *ref_E=NULL;
  v=mp->tmpF;

  for (c=0;c<C;c++){
    nt=mp->E[c]->t->size;
    assert(nt<=T);
    if (NEEDS_NORMALISATION(mp->E[c])){
      i=mp->E[c]->NormaliseByExperiment;
      ref_E=(i>=0)?(mp->E[i]):NULL;
      get_normalising_vector_with_sd(mp->E[c], ref_E, r_data, r_stdv);
      data=mp->E[c]->data;
      stdv=mp->E[c]->sd_data;      
      for (j=0;j<nt;j++){
	assert(gsl_vector_ispos(r_data[j]));
	gsl_vector_div(data[j],r_data[j]);
	//
	gsl_vector_div(stdv[j],r_data[j]);
	gsl_vector_memcpy(v,r_stdv[j]);
	assert(gsl_vector_isnonneg(v) && gsl_vector_isnonneg(data[j]));
	gsl_vector_div(v,r_data[j]);
	gsl_vector_mul(v,data[j]);
	gsl_vector_add(stdv[j],v);
      }
    }
  }
  return GSL_SUCCESS;
}
