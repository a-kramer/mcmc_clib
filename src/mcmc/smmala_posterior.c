#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include <mpi.h>
#include "smmala.h"
#include "../ode/ode_model.h"
#include "smmala_posterior.h"
#include "../app/diagnosis_output.h"
#include "model_parameters_smmala.h"
int LogLikelihood(ode_model_parameters *mp, double *l, gsl_vector *dl, gsl_matrix *FI);

int ode_solver_step(ode_solver *solver, double t, gsl_vector *y, gsl_vector* fy, gsl_matrix *yS, gsl_matrix *fyS, sensitivity_approximation *a){
  /* solves ODE for parameter vector p;
   * returns state vectors y
   */
  double tout;
  ode_model *model;
  model=solver->odeModel;
  int CVerror =  ode_solver_solve(solver, t, y->data, &tout);

  if (CVerror!=CV_SUCCESS) {
    fprintf(stderr, "ODE solver failed with ERROR = %i.\n",CVerror);
    // return a rejection message; the Likelihood is not defined for this argument;
    return GSL_EDOM;
  }
  // get sensitivities and output function values for the calculated ODE solutions
  if (ode_model_has_funcs(model)) {
    //printf("obtaining functions."); fflush(stdout);
    ode_solver_get_func(solver, tout, y->data, fy->data);
    //printf("..done\n"); fflush(stdout);
  }
  else {
    fprintf(stderr,"ode model has no output functions");
    exit(-1);
  }
  //printf("obtaining sensitivities.\n"); fflush(stdout);
  if (ode_model_has_sens(model)){
    ode_solver_get_sens(solver, tout, yS->data);
  }  else {
    // printf("approximating sensitivity\n"); fflush(stdout);
    // approximate sensitivities using steady state assumptions;
    gsl_vector_view eJt_diag;
    //gsl_permutation *permutation;
    gsl_vector_view jacp_row;
    int signum;
    int status;
    int N=ode_model_getN(model);
    int P=ode_model_getP(model);
    int i;
    // get jacobian
    ode_solver_get_jac(solver,tout,y->data,fy->data,a->jacobian_y->data);
    ode_solver_get_jacp(solver,tout,y->data,fy->data,a->jacobian_p->data);
    gsl_matrix_transpose(a->jacobian_y); // now jacobian_y(i,j)=df[i]/df[j];
    //gsl_printf("Jac_y",a->jacobian_y,1);
    gsl_matrix_memcpy(a->Jt,a->jacobian_y);
    gsl_matrix_scale(a->Jt,tout);           // this is now Jacobian*t    
    // each row of jacp is a different parameter; row1 is df[i]/dp1
    status=gsl_linalg_QR_decomp(a->jacobian_y,a->tau);
    //    status=gsl_linalg_LU_decomp(jacobian_y, permutation, &signum);
    if (status!= GSL_SUCCESS){
      fprintf(stderr,"[sens approx] QR decomposition failed. %s\n",gsl_strerror(status));      
    } 
    for (i=0;i<P;i++){
      jacp_row=gsl_matrix_row(a->jacobian_p,i);
      // solve in place; jac_p will contain the solution
      status=gsl_linalg_QR_svx(a->jacobian_y, a->tau, &jacp_row.vector);
      if (status!=GSL_SUCCESS) {
	fprintf(stderr,"[sens approx] QR solution of linear equations failed: %s. Using minimum Norm solution.\n",gsl_strerror(status));
	if (gsl_linalg_QR_lssolve(a->jacobian_y, a->tau, &jacp_row.vector, a->x, a->r)!=GSL_SUCCESS) {
	  fprintf(stderr,"[sens approx] QR lssolve also failed. exiting.\n");
	  exit(-1);
	} else {
	  gsl_vector_memcpy(&jacp_row.vector,a->x);
	}	
      }
    }
    status=gsl_linalg_exponential_ss(a->Jt,a->eJt,GSL_PREC_SINGLE);
    if (status!=GSL_SUCCESS){
      // this is not yet considered stable by GSL :/
      fprintf(stderr,"[sens approx] matrix exponential failed. %s\n",gsl_strerror(status));
    }
    eJt_diag=gsl_matrix_diagonal(a->eJt);
    gsl_vector_add_constant(&eJt_diag.vector,-1.0);
    if (gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0,a->jacobian_p,a->eJt,1.0,yS)!=GSL_SUCCESS){
      fprintf(stderr,"[sens approx] dge Matrix·Matrix failed.\n");
    }
    //gsl_permutation_free(permutation);
  } // end if has sens
  
  if (ode_model_has_funcs_sens(model)){
    ode_solver_get_func_sens(solver, tout, y->data, yS->data, fyS->data);
  } // end if has func sens

  return GSL_SUCCESS;
}


/* Normalisation information: quantities are normalised in place:
 * NORMALISATION_TYPE: fy_i ← normalisation_function(fy_i,fy_j)
 * DATA_NORMALISED_BY_REFERENCE: fy_i(t_j,k,u)←fy_i(t_j,k,u)/fy_i(t_j,k,reference_u)
 * DATA_NORMALISED_BY_TIMEPOINT: fy_i(t_j,k,u)←fy_i(t_j,k,u)/fy_i(t_l,k,u) given l
 * DATA_NORMALISED_BY_STATE_VAR: fy_i1(t_j,k,u)←fy_i1(t_j,k,u)/fy_i2(t_l,k,u) given i1 i2 and l
 * DATA_NORMALISED_BY_?:
 */

/* Note that all initial retrun values from CVODES are with respect to
 * the model's parameters k=exp(x), while we require sensitivities
 * with respect to x. Therefore, chain rule factors appear in the 
 * below calculations: dk(x)/dx = k(x); df(k(x))/d(x) = df/dk dk/dx = df/dk * k
 * 
 */
int normalise_by_reference(ode_model_parameters *mp){
  /*  normalisation by control experiment each data-set (T×N array) is
   *  normalised by the same reference set.  The reference data-set
   *  represents nominal system conditions (e.g. wild-type cells,
   *  etc.). The different data sets differ in the experimental
   *  conditions. Different conditions manifest themselves by
   *  different input-vector values:
   *            Condition_1 :     u:=[u{1,1},…,u{1,U}]
   *            Condition_2 :     u:=[u{2,1},…,u{2,U}]
   *            ...
   *            Condition_C :     u:=[u{C,1},…,u{C,U}]
   *  Reference Confition   : ref_u:=[ref_u{1},…,ref_u{U}]
   */
  int c,j,k;
  int C=mp->size->C;
  int T=mp->size->T;
  int D=mp->size->D;
  gsl_vector *fy;   
  gsl_matrix *fyS,*oS;
  gsl_vector *r_fy;
  gsl_matrix *r_fyS;
  gsl_vector *tmp;
  gsl_vector_view fyS_k_view; // 
  gsl_vector *fyS_k; // the raw fy sensitivity for one parameter k
  gsl_vector_view r_fyS_k_view; // 
  gsl_vector *r_fyS_k; // the raw reference fy sensitivity for one parameter k
  gsl_vector_view oS_k_view; // 
  gsl_vector *oS_k; // the normalised output sensitivity for one parameter k
  tmp=mp->tmpF; // maybe find a piece of unused memory somewhere?
  for (c=0;c<C;c++){
    for (j=0;j<T;j++){
      // define some shorthands
      r_fy=mp->ref_E->fy[j];
      r_fyS=mp->ref_E->fyS[j];
      fy=mp->E[c]->fy[j];
      fyS=mp->E[c]->fyS[j];
      oS=mp->E[c]->oS[j];
      // since only the first D of P parameters are sampled, only
      // those contribute to the Fisher Information, which is a
      // function of the observation sensitivity oS: the sensitivity
      // of the reproducible, normalised output, rather than the raw
      // output's fyS. So, oS will only contain derivatives with
      // respect to the first D parameters.
      gsl_vector_div(fy,r_fy);
      //printf("E[%i] t[%i]",c,j);
      //gsl_printf("fy/=r_fy",fy,0);
      
      for (k=0;k<D;k++){
	// set up vector views for the sensitivity of:
	// the raw output functions
	fyS_k_view=gsl_matrix_row(fyS,k);
	fyS_k=&(fyS_k_view.vector);
	// raw reference output functions
	r_fyS_k_view=gsl_matrix_row(r_fyS,k);
	r_fyS_k=&(r_fyS_k_view.vector);
	// normalised output (or observation)
	oS_k_view=gsl_matrix_row(oS,k);
	oS_k=&(oS_k_view.vector);
	// calculate normalised output sensitivity
	//tmp=r_fy; //r_fy is not needed anymore ! wrong
	gsl_vector_memcpy(oS_k,fyS_k);
	gsl_vector_memcpy(tmp,fy);
	gsl_vector_mul(tmp,r_fyS_k);
	gsl_vector_sub(oS_k,tmp);
	gsl_vector_div(oS_k,r_fy);
	gsl_vector_scale(oS_k,gsl_vector_get(mp->p,k));
      }
      //printf("E[%i] t[%i]",c,j);
      //gsl_printf("oS",oS,1);
    }
  }
  return GSL_SUCCESS;
}

int normalise_by_timepoint(ode_model_parameters *mp){
  /*  normalisation by the output function vector at t(l) for each
   *  experimental condition c.  Here, t(l) is the time at which the
   *  data is normalised to 1. Must be the same time for all
   *  trajectories.
   */
  int r,c,j,k,l;
  int C=mp->size->C;
  int T=mp->size->T;
  int D=mp->size->D;
  gsl_vector *fy,*l_fy;
  gsl_matrix *fyS, *oS;
  gsl_vector *tmp;
  gsl_vector_view fyS_k_view; // 
  gsl_vector *fyS_k; // the raw fy sensitivity for one parameter k
  gsl_vector_view l_fyS_k_view; // 
  gsl_vector *l_fyS_k; // the raw fy sensitivity for one parameter k
  gsl_vector_view oS_k_view; // 
  gsl_vector *oS_k; // the normalised output sensitivity for one parameter k

  tmp=mp->tmpF;
  r=mp->norm_t->size1;
  for (c=0;c<C;c++){
    l=gsl_matrix_int_get(mp->norm_t,c%r,0);    
    l_fy=mp->E[c]->fy[l];
    for (j=0;j<T;j++){
      fy=mp->E[c]->fy[j];
      fyS=mp->E[c]->fyS[j];
      oS=mp->E[c]->oS[j];
      gsl_vector_div(fy,l_fy);
      for (k=0;k<D;k++){
	// set up vector views for the sensitivity of:
	// the normalising raw output row fy at time index l
	l_fyS_k_view=gsl_matrix_row(mp->E[c]->fyS[l],k);
	l_fyS_k=&(l_fyS_k_view.vector);
	// the raw output functions
	fyS_k_view=gsl_matrix_row(fyS,k);
	fyS_k=&(fyS_k_view.vector);
	// normalised output (or observation)
	oS_k_view=gsl_matrix_row(oS,k);
	oS_k=&(oS_k_view.vector);
	// calculate normalised output sensitivity
	gsl_vector_memcpy(oS_k,fyS_k);
	gsl_vector_memcpy(tmp,fy);
	gsl_vector_mul(tmp,l_fyS_k);
	gsl_vector_sub(oS_k,tmp);
	gsl_vector_div(oS_k,l_fy);
	gsl_vector_scale(oS_k,gsl_vector_get(mp->p,k));
      }
    }	     
  }
  return GSL_SUCCESS;
}

int normalise_by_state_var(ode_model_parameters *mp){
  /*  normalisation of fy[c*T+j](i) by one of the points fy[c*T+ti](si) in the
   *  time series at t(ti), where si is a state dependent time
   *  index. So, fy[c*T+j](i)/fy[c*T+ti](si)
   */
  int c,i,j,k,ti,si;
  int rows_f, rows_t;
  int C=mp->size->C;
  int T=mp->size->T;
  int F=mp->size->F;
  int D=mp->size->D;
  gsl_vector *fy;
  gsl_vector *r_fy;
  gsl_matrix *fyS, *r_fyS, *oS;
  gsl_vector *tmp;
  gsl_vector_view si_col_view;
  gsl_vector *si_col;
  gsl_vector_view fyS_k_view; // 
  gsl_vector *fyS_k; // the raw fy sensitivity for one parameter k
  gsl_vector_view r_fyS_k_view; // 
  gsl_vector *r_fyS_k; // the raw r_fy sensitivity for one parameter k
  gsl_vector_view oS_k_view; // 
  gsl_vector *oS_k; // the normalised output sensitivity for one parameter k
  rows_f=mp->norm_f->size1;
  rows_t=mp->norm_t->size1;
  
  for (c=0;c<C;c++){
    tmp=mp->tmpF;
    r_fy=mp->ref_E->fy[0];   // these can be used as temporary storage
    r_fyS=mp->ref_E->fyS[0]; // because there is no reference experiment
    for (i=0;i<F;i++){
      si=gsl_matrix_int_get(mp->norm_f,c%rows_f,i); // normalising state index: si
      ti=gsl_matrix_int_get(mp->norm_t,c%rows_t,i); // time idx ti: fy[ti](si)
      gsl_vector_set(r_fy,i,gsl_vector_get(mp->E[c]->fy[ti],si));
      si_col_view=gsl_matrix_subcolumn(mp->E[c]->fyS[ti],si,0,D);
      si_col=&(si_col_view.vector);
      gsl_matrix_set_col(r_fyS,i,si_col);
    }
    for (j=0;j<T;j++){
      fy=mp->E[c]->fy[j];
      fyS=mp->E[c]->fyS[j];
      oS=mp->E[c]->oS[j];
      gsl_vector_div(fy,r_fy);
      for (k=0;k<D;k++){
	// set up vector views for the sensitivity of:
	// the normalising raw output row fy at time index l
	r_fyS_k_view=gsl_matrix_row(r_fyS,k);
	r_fyS_k=&(r_fyS_k_view.vector);
	// the raw output functions
	fyS_k_view=gsl_matrix_row(fyS,k);
	fyS_k=&(fyS_k_view.vector);
	// normalised output (or observation)
	oS_k_view=gsl_matrix_row(oS,k);
	oS_k=&(oS_k_view.vector);
	// calculate normalised output sensitivity
	gsl_vector_memcpy(oS_k,fyS_k);
	gsl_vector_memcpy(tmp,fy);
	gsl_vector_mul(tmp,r_fyS_k);
	gsl_vector_sub(oS_k,tmp);
	gsl_vector_div(oS_k,r_fy);
	gsl_vector_scale(oS_k,gsl_vector_get(mp->p,k));
      }
    }	     
  }
  return GSL_SUCCESS;
}


int get_normalising_vector(experiment *E, experiment *ref_E, gsl_vector *fy, int t){
  int i=E->normalised_by_experiment;
  int j,nj,k,nk;
  int f;
  double val;
  fy=NULL;
  if (ref_E!=NULL){
    assert(ref_E->normalised_by_experiment<0);
    assert(ref_E->normalised_by_output==NULL);
    assert(ref_E->normalised_by_timepoint==NULL);
  }else{
    ref_E=E;
  }
  int F=E->fy[t]->size;  
  if (E->normalised_by_timepoint!=NULL){
    nj=E->normalised_by_timepoint->size;
    if (nj==1){ // easy case: all outputs are normalised by themselves at time j
      j=gsl_vector_get(E->normalised_by_timepoint,1);
      fy=ref_E->fy[j];
    }else{ // each output is normalised differently
      fy=E->nfy;
      for (f=0;f<F;f++){
	j=gsl_vector_get(E->normalised_by_timepoint,f%nj);
	if (E->normalised_by_output!=NULL){
	  nk=ref_E->normalised_by_output->size;
	  k=gsl_vector_get(E->normalised_by_output,f%nk);
	}else{
	  k=f;
	}	 
	val=gsl_vector_get(ref_E->fy[j],k);
	gsl_vector_set(fy,f,val);
      }
    }
  }else{
    fy=ref_E->fy[t];
  }  
  return GSL_SUCCESS;
}

/* A more general normalisation routine. It reads out normalisation
 * information from the actual experiment struct and treats every
 * experiment differently, as specified. (work in progress)
 */
int normalise(ode_model_parameters *mp){
  int r,c,j,k,l;
  int C=mp->size->C;
  int T=mp->size->T;
  int D=mp->size->D;
  gsl_vector *fy,*l_fy;
  gsl_matrix *fyS, *oS;
  gsl_vector *tmp;
  gsl_vector_view fyS_k_view; // 
  gsl_vector *fyS_k; // the raw fy sensitivity for one parameter k
  gsl_vector_view l_fyS_k_view; // 
  gsl_vector *l_fyS_k; // the raw fy sensitivity for one parameter k
  gsl_vector_view oS_k_view; // 
  gsl_vector *oS_k; // the normalised output sensitivity for one parameter k
  int i;
  experiment *ref_E=NULL;
  tmp=mp->tmpF;
  r=mp->norm_t->size1;
  for (c=0;c<C;c++){
    T=mp->E[c]->t->size;
    i=mp->E[c]->normalised_by_experiment;
    ref_E=(i>0)?(mp->E[i]):NULL;
    for (j=0;j<T;j++){
      get_normalising_vector(mp->E[c], ref_E,l_fy, j);
      fy=mp->E[c]->fy[j];
      fyS=mp->E[c]->fyS[j];
      oS=mp->E[c]->oS[j];
      gsl_vector_div(fy,l_fy);
      for (k=0;k<D;k++){
	// set up vector views for the sensitivity of:
	// the normalising raw output row fy at time index l
	l_fyS_k_view=gsl_matrix_row(mp->E[c]->fyS[l],k);
	l_fyS_k=&(l_fyS_k_view.vector);
	// the raw output functions
	fyS_k_view=gsl_matrix_row(fyS,k);
	fyS_k=&(fyS_k_view.vector);
	// normalised output (or observation)
	oS_k_view=gsl_matrix_row(oS,k);
	oS_k=&(oS_k_view.vector);
	// calculate normalised output sensitivity
	gsl_vector_memcpy(oS_k,fyS_k);
	gsl_vector_memcpy(tmp,fy);
	gsl_vector_mul(tmp,l_fyS_k);
	gsl_vector_sub(oS_k,tmp);
	gsl_vector_div(oS_k,l_fy);
	gsl_vector_scale(oS_k,gsl_vector_get(mp->p,k));
      }
    }	     
  }
  return GSL_SUCCESS;
}



int assign_oS_fyS(ode_model_parameters *mp){
  int i,j,k,c;
  int C,T;
  int D,F;
  D=mp->size->D;
  F=mp->size->F;
  C=mp->size->C;
  gsl_vector *t;
  gsl_matrix_view fyS_sub;
  gsl_vector_view oS_row;
  for (c=0;c<C;c++){
    t=mp->E[c]->t;
    T=t->size;
    for (j=0;j<T;j++){
      fyS_sub=gsl_matrix_submatrix(mp->E[c]->fyS[j],0,0,D,F);
      gsl_matrix_memcpy(mp->E[c]->oS[j],&(fyS_sub.matrix));
      for (k=0;k<D;k++){
	oS_row=gsl_matrix_row(mp->E[c]->oS[j],k);
	gsl_vector_scale(&oS_row.vector,gsl_vector_get(mp->p,k));
      }
    }
  }
  return GSL_SUCCESS;
}


int LogLikelihood(ode_model_parameters *mp, double *l, gsl_vector *grad_l, gsl_matrix *fisher_information){

  /* Here, the ode integration is done
   * 
   * «l» is a scalar, the return slot of the log-likelihood
   * dl is the parameter gradient of the log likelihood (size=D)
   *
   * cvode returns functions «fy» of the state variables «y».
   * 
   * Note P is the cvode related number of ode parameters
   * while D is the number of MCMC related parameters
   * P=D+U
   */
  int i,j,c,T,C,P,U,N,D;
  double l_t_j;
  int status=GSL_SUCCESS;
  gsl_vector *y,*fy;
  gsl_matrix *yS,*fyS;
  gsl_vector *t;
  gsl_vector_view oS_row;
  gsl_vector_view input_part;
  int i_flag=GSL_SUCCESS; // ODE (i)ntegration error flag
  ode_model *model;
  ode_solver *solver;
  sensitivity_approximation *a;
  a=mp->S_approx;
  solver=mp->solver;
  model=solver->odeModel;
  C=mp->size->C;
  P=mp->size->P;
  U=mp->size->U;
  N=mp->size->N;
  D=mp->size->D;
  //int F=mp->size->F;
  //  printf("[L] D=%i\tF=%i\tU=%i\tC=%i\tP=%i\n",D,F,U,C,P);
  /* calculate reference model output if necessary:
   */
  input_part=gsl_vector_subvector(mp->p,D,U);
  //printf("[likelihood] start (normalisation type: %i).\n",mp->normalisation_type);
  //fflush(stdout);
  if (mp->normalisation_type==DATA_NORMALISED_BY_REFERENCE){
    gsl_vector_memcpy(&(input_part.vector),mp->ref_E->input_u);
    ode_solver_reinit(solver, mp->t0, mp->ref_E->init_y->data, N,
		      mp->p->data,
		      mp->p->size);
    //gsl_printf("reference p",mp->p,GSL_IS_VECTOR|GSL_IS_DOUBLE); fflush(stdout);
    //gsl_printf("ref yS0",mp->ref_E->yS0,GSL_IS_MATRIX|GSL_IS_DOUBLE);
    if (ode_model_has_sens(model)){      
      ode_solver_reinit_sens(solver, mp->ref_E->yS0->data, P, N);
    }
    t=mp->ref_E->t;
    T=t->size;
    for (j=0; j<T; j++){

      y=mp->ref_E->y[j];  //printf("y: %i, %zi\n",mp->size->N,y->size); fflush(stdout);
      fy=mp->ref_E->fy[j]; //printf("fy: %i, %zi\n",mp->size->F,fy->size);fflush(stdout);
      yS=mp->ref_E->yS[j]; //printf("yS: %i, %zi\n",mp->size->N*P,yS->size1*yS->size2);fflush(stdout);
      fyS=mp->ref_E->fyS[j]; //printf("fyS: %i, %zi\n",mp->size->F*P,fyS->size1*fyS->size2);fflush(stdout);
      i_flag&=ode_solver_step(solver, gsl_vector_get(t,j), y, fy, yS, fyS, a);
    }    
  }

  for (c=0; c<C; c++){// loop over different experimental conditions
    // write inputs into the ode parameter vector
    gsl_vector_memcpy(&(input_part.vector),mp->E[c]->input_u);
    ode_solver_reinit(solver, mp->t0, mp->E[c]->init_y->data, N,
		      mp->p->data,
		      mp->p->size);
    if (ode_model_has_sens(model)){
      ode_solver_reinit_sens(solver, mp->E[c]->yS0->data, P, N);
    }
    t=mp->E[c]->t;
    T=t->size;
    for (j=0; j<T; j++){
      y=mp->E[c]->y[j]; //printf("y: %i, %zi\n",mp->size->N,y->size);
      fy=mp->E[c]->fy[j]; //printf("fy: %i, %zi\n",mp->size->F,fy->size);
      yS=mp->E[c]->yS[j]; //printf("yS: %i, %zi\n",mp->size->N*P,yS->size1*yS->size2);
      fyS=mp->E[c]->fyS[j]; //printf("fyS: %i, %zi\n",mp->size->F*P,fyS->size1*fyS->size2);
      i_flag&=ode_solver_step(solver, gsl_vector_get(t,j), y, fy, yS, fyS, a);
    }
  }
  
  if (i_flag!=GSL_SUCCESS){
    l[0]=-INFINITY;
    return i_flag;
  }else{
    //printf("pre normalisation\n");
    //printf_omp(mp);
    
    switch (mp->normalisation_type){
    case DATA_NORMALISED_BY_TIMEPOINT:
      normalise_by_timepoint(mp);
      break;
    case DATA_NORMALISED_BY_REFERENCE:
      normalise_by_reference(mp);
      break;
    case DATA_NORMALISED_BY_STATE_VAR:
      normalise_by_state_var(mp);
      break;
    case DATA_IS_ABSOLUTE:
      assign_oS_fyS(mp);
      break;
    default:
      fprintf(stderr,"unknown normalisation method: %i\n",mp->normalisation_type);
      exit(-1);
    }
    //printf("post normalisation\n");
    //printf_omp(mp);
    
    //initialise all return values
    // log-likelihood
    l[0]=0;    
    // gradient of the log-likelihood
    gsl_vector_set_zero(grad_l);
    gsl_matrix_set_zero(fisher_information);  
  
    for (c=0; c<C; c++){
      t=mp->E[c]->t;
      T=t->size;
      for (j=0; j<T; j++){
	/* Calculate log-likelihood value
	 */
      
	gsl_vector_sub(mp->E[c]->fy[j],mp->E[c]->data[j]);
	gsl_vector_div(mp->E[c]->fy[j],mp->E[c]->sd_data[j]); // (fy-data)/sd_data
	if (gsl_blas_ddot (mp->E[c]->fy[j],mp->E[c]->fy[j], &l_t_j)!=GSL_SUCCESS){
	  printf("ddot was unsuccessful\n");
	  exit(-1);
	} // sum((fy-data)²/sd_data²)
      
	l[0]+=-0.5*l_t_j;
	/* Calculate The Likelihood Gradient and Fisher Information:
	 */
	gsl_vector_div(mp->E[c]->fy[j],mp->E[c]->sd_data[j]); // (fy-data)/sd_data²
	//gsl_dgemv(TransA,alpha,A,x,beta,y)
	status=gsl_blas_dgemv(CblasNoTrans,-1.0,mp->E[c]->oS[j],mp->E[c]->fy[j],1.0,grad_l);
	if (status!=GSL_SUCCESS){
	  gsl_printf("oS",mp->E[c]->oS[j],1);
	  gsl_printf("(fy-data)/sd_data²",mp->E[c]->fy[j],0);
	  fprintf(stderr,"dgemv was unsuccessful: %i %s\n",status,gsl_strerror(status));
	  //exit(-1);
	}
	/* Calculate the Fisher information
	 */
	for (i=0;i<D;i++) {
	  oS_row=gsl_matrix_row(mp->E[c]->oS[j],i);
	  gsl_vector_div(&(oS_row.vector),mp->E[c]->sd_data[j]);
	}
	gsl_blas_dgemm(CblasNoTrans,
		       CblasTrans, 1.0,
		       mp->E[c]->oS[j],
		       mp->E[c]->oS[j], 1.0,
		       fisher_information);	
      } // end for loop for time points
    } //end for different experimental conditions (i.e. inputs)
    //exit(0);
  }
  return i_flag & status;
}



/* Calculates the unormalised log-posterior fx, the gradient dfx, the
 * Fisher information FI.
 */
int LogPosterior(const double beta, const gsl_vector *x,  void* model_params, double *fx, gsl_vector **dfx, gsl_matrix **FI){
	
  ode_model_parameters* omp =  (ode_model_parameters*) model_params;
  double prior_value=0;            // the prior distribution related sum of squares
  int i;
  gsl_vector *prior_diff;  
  int logL_stat;
  int D=omp->prior_mu->size;

  prior_diff=omp->prior_tmp_a;

  fx[0]=0;
  for (i=0;i<D;i++) gsl_vector_set(omp->p,i,gsl_sf_exp(gsl_vector_get(x,i)));

  logL_stat=LogLikelihood(omp, &fx[1], dfx[1], FI[1]); // likelihood contribution
  //printf("likelihood: %g\n",fx[0]);
  
  gsl_vector_memcpy(prior_diff,x);
  gsl_vector_sub(prior_diff,omp->prior_mu);

  gsl_blas_dsymv(CblasUpper,-0.5,omp->prior_inverse_cov,prior_diff,0.0,dfx[2]);
  gsl_blas_ddot(prior_diff,dfx[2],&prior_value);
  fx[2]=prior_value;
  
  gsl_matrix_memcpy(FI[2],omp->prior_inverse_cov);
  int status=logL_stat;
  //printf("PosteriorFI: prior_value=%f\n",prior_value);

  fx[0]=fx[1];
  fx[0]*=beta;
  fx[0]+=fx[2];

  gsl_vector_memcpy(dfx[0],dfx[1]);
  gsl_vector_scale(dfx[0],beta);
  gsl_vector_add(dfx[0],dfx[2]);

  gsl_matrix_memcpy(FI[0],FI[1]);
  gsl_matrix_scale(FI[0],beta*beta);
  gsl_matrix_add(FI[0],FI[2]);
  //printf("[LogPosterior] done.\n"); fflush(stdout);
  return status;
}


