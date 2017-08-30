/*
 *  ode_smmala.c
 *
 *  Copyright (C) 2012,2013,2015,2016 Vassilios Stathopoulos stathv@gmail.com, Andrei
 *  Kramer andrei.kramer@scilifelab.se
 *
 *	This program is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 */


/* Nomenclature: some variables are not self-descriptive. Here we
 * mostly adhere to CVODE's variable names: y is the numerical
 * solution to an initial value problem (IVP). The IVP has unknown
 * parameters k=exp(x).  Here, x is the MCMC sampling
 * variable. Sampling is done in log-space. So, exp(x) is passed to
 * the IVP as parameters. The model has internal "Expressions" which
 * remain unnamed here and output functions (measurable quantites),
 * which are named fy. Both y and fy have sensitivities with respect
 * to the model parameters (and sampling parameters) yS=dy(t)/dk and
 * similarly for fyS. Some variables are converted "in place", so
 * their interpretation might change slightly during intermediate
 * steps.
 *
 * In the realm of systems biology, the typical variable names are: 
 * x_i - state variables of the ODE model
 * k_j - parameters of the model
 * These terminologies might be mixed up in the comments [TODO: clean up comments]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <mpi.h>
#include "read_cnf.h"
#include "../mcmc/smmala.h"
#include "../ode/ode_model.h"


// sampling actions 
#define SMPL_RESUME 1
#define SMPL_FRESH 0
#define SMPL_RESUME_TUNE 2


/* Auxiliary structure with working storage and aditional parameters for
 * a multivariate normal model with known covariance matrix and zero mean.
 *
typedef struct {
	int D;
	double* Variance;
	double* Precision;
	double* tmpVec;
	char init;
} mvNormParams;
 */
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


/* Calculates the unormalised posterior fx, the gradient dfx, the Fisher information FI 
 * and the partial derivatives of the Fisher information matrix for a multivariate normal.
 * Arguments:
 *  x:				values for the model parameters.
 *	model_params:	a pointer to mvNormParams struct.
 * Returns:
 *  fx:				the unormalised posterior density at x.
 *  dfx:			the gradient of fx w.r.t x. A double array with size N equal to the size of x.
 *  FI:				the Fisher information evaluated at x. A double array with size N*N whith 
 *					elements stored in row major order. The i,j element of the FI is stored at FI[i*N + j].
 *  dFI:			the partial derivatives of FI w.r.t. x. A double array with size N*N*N.
 *					The i,j element of the k^th partial derivative of FI, \frac{\partial FI}{\partial x_k}
 *					is stored at dFI[k*N*N + i*N + j].
 * The function should return 0 for success or a non-zero value indicating failure.
 */
int Posterior(const double* x,  void* model_params, double* fx, double* dfx, double* FI);
//int PriorRnd(gsl_rng* rng, const const void* model_params, double* sample);

int ode_solver_step(ode_solver *solver, double t, gsl_vector *y, gsl_vector* fy, gsl_matrix *yS, gsl_matrix *fyS){
  /* solves ODE for parameter vector p;
   * returns state vectors y
   */
  double tout;
  ode_model *model;
  model=solver->odeModel;
  int CVerror =  ode_solver_solve(solver, t, y->data, &tout);

  if (CVerror) {
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
    gsl_matrix *jacobian_y;
    gsl_matrix *jacobian_p;
    gsl_matrix *eJt;
    gsl_matrix *Jt;
    gsl_vector *tau;
    gsl_vector *x;
    gsl_vector *r;
    gsl_vector_view eJt_diag;
    //gsl_permutation *permutation;
    gsl_vector_view jacp_row;
    int signum;
    int status;
    int N=ode_model_getN(model);
    int P=ode_model_getP(model);
    int i;
    //permutation=gsl_permutation_alloc(N);
    jacobian_y=gsl_matrix_alloc(N,N);
    jacobian_p=gsl_matrix_alloc(P,N);
    tau=gsl_vector_alloc(N);
    x=gsl_vector_alloc(N);
    r=gsl_vector_alloc(N);
    eJt=gsl_matrix_alloc(N,N);
    Jt=gsl_matrix_alloc(N,N);
    
    // get jacobian
    ode_solver_get_jac(solver,tout,y->data,fy->data,jacobian_y->data);
    ode_solver_get_jacp(solver,tout,y->data,fy->data,jacobian_p->data);
    gsl_matrix_transpose(jacobian_y); // now jacobian_y(i,j)=df[i]/df[j];
    //gsl_printf("Jac_y",jacobian_y,1);
    gsl_matrix_memcpy(Jt,jacobian_y);
    gsl_matrix_scale(Jt,tout);           // this is now Jacobian*t    
    // each row of jacp is a different parameter; row1 is df[i]/dp1
    status=gsl_linalg_QR_decomp(jacobian_y,tau);
    //    status=gsl_linalg_LU_decomp(jacobian_y, permutation, &signum);
    if (status!= GSL_SUCCESS){
      fprintf(stderr,"[sens approx] QR decomposition failed. %s\n",gsl_strerror(status));      
    } 
    for (i=0;i<P;i++){
      jacp_row=gsl_matrix_row(jacobian_p,i);
      // solve in place; jac_p will contain the solution
      status=gsl_linalg_QR_svx(jacobian_y, tau, &jacp_row.vector);
      if (status!=GSL_SUCCESS) {
	fprintf(stderr,"[sens approx] QR solution of linear equations failed: %s. Using minimum Norm solution.\n",gsl_strerror(status));
	if (gsl_linalg_QR_lssolve(jacobian_y, tau, &jacp_row.vector, x, r)!=GSL_SUCCESS) {
	  fprintf(stderr,"[sens approx] QR lssolve also failed. exiting.\n");
	  exit(-1);
	} else {
	  gsl_vector_memcpy(&jacp_row.vector,x);
	}	
      }
    }
    status=gsl_linalg_exponential_ss(Jt,eJt,GSL_PREC_SINGLE);
    if (status!=GSL_SUCCESS){
      // this is not yet considered stable by GSL :/
      fprintf(stderr,"[sens approx] matrix exponential failed. %s\n",gsl_strerror(status));
    }
    eJt_diag=gsl_matrix_diagonal(eJt);
    gsl_vector_add_constant(&eJt_diag.vector,-1.0);
    if (gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0,jacobian_p,eJt,1.0,yS)!=GSL_SUCCESS){
      fprintf(stderr,"[sens approx] dge Matrix·Matrix failed.\n");
    }
    //gsl_permutation_free(permutation);
    gsl_matrix_free(jacobian_y);
    gsl_matrix_free(jacobian_p);
    gsl_vector_free(tau);
    gsl_vector_free(x);
    gsl_vector_free(r);
    gsl_matrix_free(Jt);
    gsl_matrix_free(eJt);
  } // end if has sens
  
  if (ode_model_has_funcs_sens(model)){
    ode_solver_get_func_sens(solver, tout, y->data, yS->data, fyS->data);
  } // end if has func sens

  return GSL_SUCCESS;
}

int printf_omp(ode_model_parameters *mp){
  int i,j,k,c,T,D,C;
  gsl_vector *t;
  char name[512];

  C=mp->size->C;
  D=mp->size->D;
  gsl_printf("p",mp->p,0);
  printf("Reference Experiment\n");
  gsl_printf("ref y init",mp->ref_E->init_y,0);
  t=mp->ref_E->t;
  if (t!=NULL){
    T=t->size;
    for (j=0; j<T; j++){
      printf("t[%i]\n",j);
      sprintf(name,"ref fy(t[%i])",j);
      //gsl_printf(sprintf("ref y(%i)",j),mp->ref_E->y[j],0);
      gsl_printf(name,mp->ref_E->fy[j],0);
      //gsl_printf(sprintf("ref yS(%i)",j),mp->ref_E->yS[j],1);
      sprintf(name,"ref fyS(t[%i])",j);
      gsl_printf(name,mp->ref_E->fyS[j],1);
    }
  }
  for (c=0;c<C;c++){
    printf("Experiment %i\n",c);
    gsl_printf("y_init",mp->E[c]->init_y,0);
    t=mp->E[c]->t;
    if (t!=NULL){
    T=t->size;
    for (j=0; j<T; j++){
      printf("t[%i]\n",j);
      //gsl_printf("y",mp->E[c]->y[j],0);
      sprintf(name,"E[%i] fy(t[%i])",c,j);
      gsl_printf(name,mp->E[c]->fy[j],0);
      sprintf(name,"E[%i] data(t[%i])",c,j);
      gsl_printf(name,mp->E[c]->data[j],0);
      sprintf(name,"E[%i] sd data(t[%i])",c,j);
      gsl_printf(name,mp->E[c]->sd_data[j],0);
      //gsl_printf("yS",mp->E[c]->yS[j],1);
      sprintf(name,"E[%i] fyS(t[%i])",c,j);
      gsl_printf(name,mp->E[c]->fyS[j],1);
      sprintf(name,"E[%i] oS(t[%i])",c,j);
      gsl_printf(name,mp->E[c]->oS[j],1);
    }
    } else {
      printf("E[%i] not properly defined; t missing\n",c);
    }    
  }
  fflush(stdout);
  return EXIT_SUCCESS;
}


int LikelihoodComplexNorm(ode_model_parameters *mp, double *l, double *dl, double *FI){

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
  int status;
  gsl_vector *y,*fy;
  gsl_matrix *yS,*fyS;
  gsl_vector *t;
  gsl_vector_view grad_l_view;
  gsl_vector *grad_l;
  gsl_matrix_view fisher_information_view;
  gsl_matrix *fisher_information;
  gsl_vector_view oS_row;
  gsl_vector_view input_part;
  int i_flag=CV_SUCCESS;
  ode_model *model;
  ode_solver *solver;

  solver=mp->solver;
  model=solver->odeModel;
  C=mp->size->C;
  P=mp->size->P;
  U=mp->size->U;
  N=mp->size->N;
  D=mp->size->D;
  //  printf("[L] D=%i\tF=%i\tU=%i\tC=%i\tT=%i\tP=%i\n",D,F,U,C,T,P);
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
    //gsl_printf("reference p",mp->p,0); fflush(stdout);
    //gsl_printf("ref yS0",mp->ref_E->yS0,1);
    if (ode_model_has_sens(model)){      
      ode_solver_reinit_sens(solver, mp->ref_E->yS0->data, P, N);
    }
    t=mp->ref_E->t;
    T=t->size;
    for (j=0; j<T; j++){

      y=mp->ref_E->y[j];  //printf("y: %i, %i\n",mp->size->N,y->size); fflush(stdout);
      fy=mp->ref_E->fy[j]; //printf("fy: %i, %i\n",mp->size->F,fy->size);fflush(stdout);
      yS=mp->ref_E->yS[j]; //printf("yS: %i, %i\n",mp->size->N*P,yS->size1*yS->size2);fflush(stdout);
      fyS=mp->ref_E->fyS[j]; //printf("fyS: %i, %i\n",mp->size->F*P,fyS->size1*fyS->size2);fflush(stdout);
      ode_solver_step(solver, gsl_vector_get(t,j), y, fy, yS, fyS);
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
      y=mp->E[c]->y[j]; //printf("y: %i, %i\n",mp->size->N,y->size);
      fy=mp->E[c]->fy[j]; //printf("fy: %i, %i\n",mp->size->F,fy->size);
      yS=mp->E[c]->yS[j]; //printf("yS: %i, %i\n",mp->size->N*P,yS->size1*yS->size2);
      fyS=mp->E[c]->fyS[j]; //printf("fyS: %i, %i\n",mp->size->F*P,fyS->size1*fyS->size2);
      i_flag=ode_solver_step(solver, gsl_vector_get(t,j), y, fy, yS, fyS);
    }
  }
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
  grad_l_view=gsl_vector_view_array(dl,D);
  grad_l=&(grad_l_view.vector);
  gsl_vector_set_zero(grad_l);
  // fisher information
  fisher_information_view=gsl_matrix_view_array(FI,D,D);
  fisher_information=&(fisher_information_view.matrix);
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
	printf("dgemv was unsuccessful: %i %s\n",status,gsl_strerror(status));
	exit(-1);
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
  return i_flag;
}
int print_help(){
  printf("Usage:\n");
  printf("-c ./data.cfg\n");
  printf("\t\t\tdata.cfg contains the data points and the conditions of measurement.\n\n");
  printf("-l ./ode_model.so\n");
  printf("\t\t\tode_model.so is a shared library containing the CVODE functions of the model.\n\n");
  printf("-s $N\n");
  printf("\t\t\t$N sample size. default N=10.\n\n");
  printf("-r, --resume\n");
  printf("\t\t\tresume from last sampled MCMC point. Only the last MCMC position is read from the file named «resume.double». Everything else about the problem can be changed.\n");
  printf("-o ./output_file\n");
  printf("\t\t\tfile for output. Output can be binary.\n\n");
  printf("-b\n");
  printf("\t\t\toutput mode: binary.\n\n");
  printf("-a $a\n");
  printf("\t\t\ttarget acceptance value (markov chain will be tuned for this acceptance).\n\n");
  printf("--seed $seed\n");
  printf("\t\tset the gsl pseudo random number generator seed to $seed.\n\n");
  
  //  printf("test for 1/Inf=%f\n",1.0/INFINITY);
  return EXIT_SUCCESS;
}

void print_chunk_graph(gsl_matrix *X, gsl_vector *lP){
  int width=100; // we assume that the display can show $width characters
  int i,j,k,n,nc;
  int tmp;
  double *x;
  gsl_vector_view x_view;
  double Q[5];
  int q[5];
  double max,min,range;
  char s[32],c[32];
  n=X->size2;
  max=gsl_matrix_max(X);
  min=gsl_matrix_min(X);
  range=max-min;
  printf("range: [%g,%g] (%g)\n",min,max,range);
  for (i=0;i<X->size1;i++){
    //sort each row:
    x_view=gsl_matrix_row(X,i);
    gsl_sort_vector(&(x_view.vector));
    //determine eachquantile
    x=gsl_matrix_ptr(X,i,0);
    Q[0]=gsl_stats_quantile_from_sorted_data(x,1,n,0.01);
    Q[1]=gsl_stats_quantile_from_sorted_data(x,1,n,0.25);
    Q[2]=gsl_stats_quantile_from_sorted_data(x,1,n,0.50);
    Q[3]=gsl_stats_quantile_from_sorted_data(x,1,n,0.75);
    Q[4]=gsl_stats_quantile_from_sorted_data(x,1,n,0.99);
    //printf("quantiles: %g\t%g\t%g\t%g\t%g\n",Q[0],Q[1],Q[2],Q[3],Q[4]);

    for (j=0;j<5;j++) {
      q[j]=(int) ((Q[j]-min)*width/range);
    }
    sprintf(s," -LU- ");
    sprintf(c,"+{|}+ ");
    tmp=0;
    for (k=0;k<5;k++){
      nc=q[k]-tmp;
      for (j=0;j<nc;j++) {
	printf("%c",s[k]);
      }
      tmp=q[k];
      printf("%c",c[k]);
    }
    printf("\n\n");    
  }
  printf("|");
  for (j=0;j<width-2;j++) printf("-"); printf("|\n");
  printf("%+4.4g",min);
  for (j=0;j<width-8;j++) printf(" ");
  printf("%+4.4g\n",max);  
}



int main (int argc, char* argv[]) {
  int D = 0; // number of MCMC sampling variables, i.e. model parameters
  //int C = 0; // number of experimental conditions, i.e different input vectors
  int i=0;
  int warm_up=0; // sets the number of burn in points at command line
  char *cfilename=NULL;
  char lib_name[512];
  ode_model_parameters omp;
  FILE *cnf;   // configuration file, with file name: cfilename
  char global_sample_filename_stem[512]="sample.dat"; // filename basis
  char rank_sample_file[512]; // filename for sample output
  //char *x_sample_file=NULL; // filename for sample output x(t,p)
  //char *y_sample_file=NULL; // filename for sample output y(t,p)
  FILE *oFile; // will be the file named «sample_file»
  FILE *rFile; // last sampled value will be written to this file
  char resume_filename[512]="resume.double";
  int output_is_binary=0;
  double seed = 1;
  int sampling_action=SMPL_FRESH;
  size_t resume_count;
  gsl_error_handler_t *gsl_error_handler;
  int start_from_prior=0;
  gsl_error_handler = gsl_set_error_handler_off();
  int sensitivity_approximation=0;
  //  strcpy(sample_file,"sample.dat");
  main_options cnf_options;
  cnf_options.output_file=global_sample_filename_stem;
  cnf_options.library_file=lib_name;
  cnf_options.target_acceptance=-0.5;
  cnf_options.initial_stepsize =-0.1; // not used here in smmala
  cnf_options.sample_size=-10;
  MPI_Init(&argc,&argv);
  int rank,R;
  MPI_Comm_size(MPI_COMM_WORLD,&R);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  for (i=0;i<argc;i++){
    if (strcmp(argv[i],"-c")==0) cfilename=argv[i+1];
    else if (strcmp(argv[i],"-p")==0 || strcmp(argv[i],"--prior-start")==0) start_from_prior=1;
    else if (strcmp(argv[i],"-w")==0 || strcmp(argv[i],"--warm-up")==0) warm_up=strtol(argv[i+1],NULL,10);
    else if (strcmp(argv[i],"--resume")==0 || strcmp(argv[i],"-r")==0) sampling_action=SMPL_RESUME;
    else if (strcmp(argv[i],"--sens-approx")==0) sensitivity_approximation=1;
    else if (strcmp(argv[i],"-l")==0) strcpy(cnf_options.library_file,argv[i+1]);
    //    else if (strcmp(argv[i],"-n")==0) Tuning=0;
    else if (strcmp(argv[i],"-s")==0) cnf_options.sample_size=strtol(argv[i+1],NULL,0);
    else if (strcmp(argv[i],"-o")==0) strcpy(cnf_options.output_file,argv[i+1]);
    else if (strcmp(argv[i],"-b")==0) output_is_binary=1;
    else if (strcmp(argv[i],"-a")==0) cnf_options.target_acceptance=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-i")==0) cnf_options.initial_stepsize=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"--seed")==0) seed=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-h")==0) {print_help(); exit(0);}
    
    //else printf("unknown option {%s}.\n",argv[i]);
  }

  seed=seed*137+13*rank;

  /* load model */
  ode_model *odeModel = ode_model_loadFromFile(lib_name);  /* alloc */
  if (odeModel == NULL) {
    fprintf(stderr, "# Library %s could not be loaded.\n",lib_name);
    exit(1);
  } else printf( "# Library %s loaded.\n",lib_name);

  sprintf(resume_filename,"%s_resume_%04l.double",lib_name,rank);
  if (strcmp(cnf_options.output_file,global_sample_filename_stem)==0)
    sprintf(rank_sample_file,"rank_%04l_of_%l_%s_%s",rank,R,lib_name,global_sample_filename_stem);
  else
    sprintf(rank_sample_file,"rank_%04l_of_%l_%s",rank,R,cnf_options.output_file);
  cnf_options.output_file=rank_sample_file;

  ode_solver* solver = ode_solver_alloc(odeModel); /* alloc */
  if (solver == NULL) {
    fprintf(stderr, "# Solver %s could not be created.\n",lib_name);
    ode_model_free(odeModel);
    exit(1);
  } else printf("# Solver for %s created.\n",lib_name);

  if (sensitivity_approximation){
    ode_solver_disable_sens(solver);
    odeModel->vf_sens=NULL; // make sens function unavailable; ode_model_has_sens(model) will return «False»;
  }
  /* init solver */
  realtype solver_param[3] = {ODE_SOLVER_ABS_ERR, ODE_SOLVER_REL_ERR, 0};
  
  /* define local variables for parameters and inital conditions */
  int N = ode_model_getN(odeModel);
  int P = ode_model_getP(odeModel);
  int F = ode_model_getF(odeModel);
  omp.size=(problem_size*) malloc(sizeof(problem_size));
  omp.size->N=N;
  omp.size->P=P;
  omp.size->F=F;

  double y[N];
  gsl_vector_view y_view=gsl_vector_view_array(y,N);
	
  double p[P];
  gsl_vector_view p_view=gsl_vector_view_array(p,P);
  ode_model_get_default_params(odeModel, p, P);
  //gsl_printf("default parameters",&(p_view.vector),0);
  //  double t_old = solver_param[2];

  omp.solver=solver;

  cnf=fopen(cfilename,"r");
  if (cnf!=NULL){
    printf("# reading configuration.\n");
    parse_config(cnf,&omp,&cnf_options);
  } else {
    fprintf(stderr,"# Could not open config file %s.\n",cfilename);
    exit(1);
  }
  fclose(cnf);
 
  cnf_options.initial_stepsize=fabs(cnf_options.initial_stepsize);
  cnf_options.target_acceptance=fabs(cnf_options.target_acceptance);
  cnf_options.sample_size=fabs(cnf_options.sample_size);

  if (omp.normalisation_type==DATA_IS_ABSOLUTE){
    ode_model_get_initial_conditions(odeModel, y, N);
  } else {
    gsl_vector_set_all(&(y_view.vector),1.0);
  }

  
  /* gsl_printf("data",omp.Data,1); */
  /* gsl_printf("sd data",omp.sdData,1); */
  /* gsl_printf("u",omp.input_u,1); */
  /* gsl_printf("t",omp.t,0); */
 
  // unspecified initial conditions
  if (omp.ref_E->init_y==NULL){
    omp.ref_E->init_y=gsl_vector_alloc(N);
    gsl_vector_memcpy(omp.ref_E->init_y,&(y_view.vector));
  }
  
  for (i=0;i<omp.size->C;i++){
    if (omp.E[i]->init_y==NULL){
      omp.E[i]->init_y=gsl_vector_alloc(N);
      gsl_vector_memcpy(omp.E[i]->init_y,&(y_view.vector));
    }
  }  
  printf("# init ivp: t0=%g\n",omp.t0); fflush(stdout);
  ode_solver_init(solver, omp.t0, y, N, p, P);
  printf("# solver initialised.\n"); fflush(stdout);
  ode_solver_setErrTol(solver, solver_param[1], &solver_param[0], 1);
  if (ode_model_has_sens(odeModel)) {
    ode_solver_init_sens(solver, omp.ref_E->yS0->data, P, N);
    printf("# sensitivity analysis initiated.\n");
  }

  printf("# allocating memory for the SMMALA model.\n");
  fflush(stdout);
  smmala_model* model = smmala_model_alloc(Posterior, NULL, &omp); // ode_model_parameters
  /* initial parameter values */
  D=omp.size->D;
  double init_x[D];
  /* allocate a new RMHMC MCMC kernel */
  printf("# allocating memory for a new SMMALA MCMC kernel.\n");
  mcmc_kernel* kernel = smmala_kernel_alloc(D,
					    cnf_options.initial_stepsize,
					    model,
					    seed,
					    cnf_options.target_acceptance);
  /* initialise MCMC */
  //gsl_printf("prior mean",omp.prior_mu,0);  
  if (sampling_action==SMPL_RESUME){
    rFile=fopen(resume_filename,"r");
    if (rFile==NULL) {
      printf("Could not open resume file. Starting from: ");
      for (i=0;i<D;i++) init_x[i]=gsl_vector_get(omp.prior_mu,i);
    } else {
      resume_count=fread(init_x, sizeof(double), D, rFile);
      fclose(rFile);
      if (resume_count!=D) {
	fprintf(stderr,"Reading from resume file returned a wrong number of values.");
	exit(-1);
      }
    }
  } else if (start_from_prior==1){
    printf("# setting initial mcmc vector to prior mean.\n");
    for (i=0;i<D;i++) init_x[i]=gsl_vector_get(omp.prior_mu,i);
  } else {
    printf("# setting mcmc initial value to log(default parameters)\n");
    for (i=0;i<D;i++) init_x[i]=gsl_sf_log(p[i]);
  }
  
  printf("# initializing MCMC.\n");
  printf("init_x:");
  for (i=0;i<D;i++) printf(" %g ",init_x[i]);
  printf("\n");
  mcmc_init(kernel, init_x);

  //gsl_printf("Normalised Data",omp.Data,1);
  //gsl_printf("standard deviation",omp.sdData,1); //exit(0);
  //fflush(stdout);

  printf("# test evaluation of Posterior function.\n");
  fflush(stdout);
  /*inits*/
  double fx;
  double dfx[D];
  double FI[D*D];
  Posterior(init_x, &omp, &fx, dfx, FI);
  printf("# θ=θ₀; Posterior(θ|D)=%g;\n",fx);

  /* print first sample, initial values in init_x */
  mcmc_print_sample(kernel, stdout);
  ode_solver_print_stats(solver, stdout);
  size_t acc_c = 0;
  double acc_rate;
  size_t it;
  int acc;
	
  oFile=fopen(cnf_options.output_file,"w");
  if (oFile!=NULL){
    printf("# writing sample to output file.\n");
  } else {
    fprintf(stderr,"# Could not open output file {%s} for writing.\n",cnf_options.output_file);
    exit(1);
  }

  size_t BurnInSamples;
  if (warm_up==0){
    BurnInSamples = 7 * (int) sqrt(cnf_options.sample_size);
  } else {
    BurnInSamples=warm_up;
  }

  if (sampling_action==SMPL_FRESH){
    fprintf(stdout, "# Burn-in phase.\n");
    gsl_vector *lP;
    gsl_matrix *X;
    gsl_vector_view X_view, array_view;
    size_t Chunk=100;
    int beta_bar_width=60; //characters, on screen
    int b,B;
    lP=gsl_vector_alloc(Chunk);
    X=gsl_matrix_alloc(D,Chunk);
  
    /* "find mode" Loop */
    for (it = 0; it < BurnInSamples; it++) {
      omp.beta=((double) it)/((double) BurnInSamples);
      omp.beta=gsl_pow_6(omp.beta);
      omp.beta=omp.beta/(gsl_pow_6(0.5)+omp.beta);
      mcmc_sample(kernel, &acc);
      acc_c += acc;
      /* print sample */
      //mcmc_print_sample(kernel, stdout);
      //ptr=gsl_matrix_const_ptr(X,i,0);
      X_view=gsl_matrix_column(X,it);
      array_view=gsl_vector_view_array(kernel->x,D);
      gsl_vector_memcpy(&(X_view.vector),&(array_view.vector));
      //printf("accessing fx: %g.\n",fx);
      fx=kernel->fx[0];
      //printf("accessing fx: %g.\n",fx);
      gsl_vector_set(lP,it,fx);

      /* Addapt MCMC parameters every 100 burn-in samples */
      if ( ((it + 1) % Chunk) == 0 ) {
	acc_rate = (double) acc_c / (double) Chunk;
	fprintf(stdout, "# Iteration: %li\tAcceptance rate: %.2g\t",it, acc_rate);
	/* print log and statistics */
	mcmc_print_stats(kernel, stdout);
	printf("beta: [\033[22;32m");
	B=beta_bar_width;
	b=B*omp.beta;
	for (i=0;i<b;i++) printf("#");
	for (i=b;i<B;i++) printf("¤");
	printf("\033[0m]\n");
	print_chunk_graph(X,lP);
	/* adapt parameters */
	mcmc_adapt(kernel, acc_rate);
	acc_c = 0;
	fflush(stdout);
      }
    }

    omp.beta=1.0;

    /* Burn In Loop and Tuning*/
    for (it = 0; it < BurnInSamples; it++) {
      mcmc_sample(kernel, &acc);
      acc_c += acc;
      /* print sample */
      //mcmc_print_sample(kernel, stdout);
      if ( ((it + 1) % 100) == 0 ) {
	acc_rate = (double)acc_c / (double)100;
	fprintf(stdout, "# Iteration: %li\tAcceptance rate: %.2g\t",it, acc_rate);
	/* print log and statistics */
	mcmc_print_stats(kernel, stdout);
	/* adapt  parameters */
	mcmc_adapt(kernel, acc_rate);
	acc_c = 0;
      }
    }
    fprintf(stdout, "\n# Burn-in complete, sampling from the posterior.\n");
  }

  /* full Posterior loop */
  omp.beta=1.0;
  size_t Samples = cnf_options.sample_size;
  clock_t ct=clock();
  
  for (it = 0; it < Samples; it++) {
    /* draw a sample using RMHMC */
    mcmc_sample(kernel, &acc);
    acc_c += acc;
		
    /* print sample */
    if (output_is_binary) {
      mcmc_write_sample(kernel, oFile);
    }else{
      mcmc_print_sample(kernel, oFile);
    }
    /* print sample log and statistics every 100 samples */
    if ( ((it + 1) % 100) == 0 ) {
      acc_rate = (double)acc_c / (double)100;
      fprintf(stdout, "# Iteration: %li\tAcceptance rate: %.2g\t",it, acc_rate);
      mcmc_print_stats(kernel, stdout);
      acc_c = 0;
    }
  }
  rFile=fopen(resume_filename,"w");
  mcmc_write_sample(kernel, rFile);
  fclose(rFile);
  ct=clock()-ct;
  fclose(oFile);
  printf("# computation time spend sampling: %f s\n",((double) ct)/((double) CLOCKS_PER_SEC));


  /* clear memory */
  smmala_model_free(model);
  mcmc_free(kernel);
  ode_model_parameters_free(&omp);
  MPI_Finalize();
  return EXIT_SUCCESS;
}



/* Calculates the unormalised log-posterior fx, the gradient dfx, the
 * Fisher information FI.
 */
int Posterior(const double* x,  void* model_params, double* fx, double* dfx, double* FI){
	
  ode_model_parameters* omp =  (ode_model_parameters*) model_params;
  double prior_value=0;            // the prior distribution related sum of squares
  int i;
  gsl_matrix *inv_cov;
  gsl_vector *prior_diff, *prior_mv;
  
  const gsl_vector *x_v;
  // double l;
  int logL_stat;
  int D=omp->size->D;

  prior_diff=omp->prior_tmp_a;
  prior_mv=omp->prior_tmp_b;
  inv_cov=omp->prior_inverse_cov;
  
  gsl_vector_const_view x_view=gsl_vector_const_view_array(x,D);
  x_v=&x_view.vector;
  
  fx[0]=0;
  for (i=0;i<D;i++) gsl_vector_set(omp->p,i,gsl_sf_exp(x[i]));
  logL_stat=LikelihoodComplexNorm(omp, fx, dfx, FI);
  //printf("likelihood: %g\n",fx[0]);

  gsl_vector_memcpy(prior_diff,x_v);
  gsl_vector_sub(prior_diff,omp->prior_mu);

  gsl_blas_dsymv(CblasUpper,-0.5,inv_cov,prior_diff,0.0,prior_mv);
  gsl_blas_ddot(prior_diff,prior_mv,&prior_value);
  //  printf("PosteriorFI: prior_value=%f\n",prior_value);
  fx[0]*=omp->beta;
  fx[0]+=prior_value;
  /* here the "sum_of_squares" belongs to the log-likelihood
   * while the prior_value is the prior-related part of the log-posterior
   * these values omit all constants.
   */

  // now for the gradient:
  for (i=0; i < D ; i++){
    dfx[i]*= omp->beta;
    dfx[i]+= gsl_vector_get(prior_mv,i); // prior component
  }
  
  /* FI is constant, it does not depend on the parameters x ... */
  for (i=0; i < D*D ; i++){
    FI[i] *= gsl_pow_2(omp->beta);
    FI[i] += omp->prior_inverse_cov->data[i];
  }
  return logL_stat;
}
