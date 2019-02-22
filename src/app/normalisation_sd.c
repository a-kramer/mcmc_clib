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
int get_normalising_vector_with_sd(experiment *E, experiment *ref_E){
  //int i=E->NormaliseByExperiment;
  int t,rt;
  int f,k,nk;
  double val, dval;
  int j=E->NormaliseByTimePoint;
  
  if(ref_E==NULL){
    ref_E=E;
    assert(j>=0);
  }
  int F=E->data[0]->size;
  int T=E->t->size;
  int rT=ref_E->t->size;

  assert(j<T);
  if (j<0) assert(T==rT);  
  if (E->NormaliseByOutput==NULL){  
    for (t=0;t<T;t++){
      // find the right reference data and stdv
      rt=j<0?t:j;
      gsl_vector_memcpy(E->normalise->data[t],ref_E->data[rt]);
      gsl_vector_memcpy(E->normalise->stdv[t],ref_E->sd_data[rt]);      
    }
  }else{ // each output is normalised differently
    // copy elements to fy and fyS
    for (f=0;f<F;f++){
      // output function:
      k=gsl_vector_int_get(E->NormaliseByOutput,f)%F;
      for (t=0;t<T;t++){
	rt=j<0?t:j;      
	val=gsl_vector_get(ref_E->data[rt],k);
	dval=gsl_vector_get(ref_E->sd_data[rt],k);
	gsl_vector_set(E->normalise->data[t],f,val);
	gsl_vector_set(E->normalise->stdv[t],f,dval);
      }
    }
  }  
  return GSL_SUCCESS;
}

int normalise_with_sd(void *model_parameters){
  ode_model_parameters *mp=model_parameters;
  int c,t,k,l;
  int C=mp->size->C;
  int T=mp->size->T;
  int D=mp->size->D;
  int F=mp->size->F;
  assert(T>0);
  gsl_vector **data, **stdv;
  gsl_vector **r_data, **r_stdv;
  gsl_vector *v;
  int i,nt;
  experiment *ref_E=NULL;
  v=mp->tmpF;

  for (c=0;c<C;c++){
    nt=mp->E[c]->t->size;
    assert(nt<=T);
    if (NEEDS_NORMALISATION(mp->E[c])){
      i=mp->E[c]->NormaliseByExperiment;
      ref_E=(i<0)?NULL:(mp->E[i]);
      get_normalising_vector_with_sd(mp->E[c], ref_E);
      data=mp->E[c]->data;
      stdv=mp->E[c]->sd_data;
      r_data=mp->E[c]->normalise->data;
      r_stdv=mp->E[c]->normalise->stdv;
      for (t=0;t<nt;t++){
	assert(gsl_vector_ispos(r_data[t]));
	gsl_vector_div(data[t],r_data[t]);
	//
	gsl_vector_div(stdv[t],r_data[t]);
	gsl_vector_memcpy(v,r_stdv[t]);
	assert(gsl_vector_isnonneg(v) && gsl_vector_isnonneg(data[t]));
	gsl_vector_div(v,r_data[t]);
	gsl_vector_mul(v,data[t]);
	gsl_vector_add(stdv[t],v);
      }
    }
  }
  return GSL_SUCCESS;
}
