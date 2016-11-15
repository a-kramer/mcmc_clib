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
#include "read_cnf.h"
#include "../mcmc/smmala.h"
#include "../ode/ode_model.h"

// sampling actions 
#define SMPL_RESUME 1
#define SMPL_FRESH 0
#define SMPL_RESUME_TUNE 2

// normalisation methods
#define DATA_NORMALISED_BY_REFERENCE 0
#define DATA_NORMALISED_BY_TIMEPOINT 1
#define DATA_NORMALISED_BY_STATE_VAR 2


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
  
  printf("[%s]",name);
  if (is_matrix){
    A=(gsl_matrix*) gsl_thing;
    printf(" %i × %i elements\n",A->size1,A->size2);
    gsl_matrix_fprintf(stdout, A, "%g");
  }else{
    x=(gsl_vector*) gsl_thing;
    printf(" %i elements\n",x->size);
    gsl_vector_fprintf(stdout, x, "%g");
  }
  printf("[/%s]\n",name);
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
  int c,i,j,k;
  int C=mp->input_u->size1;
  int T=mp->t->size1;
  int F=mp->y[0]->size2;
  int D=mp->D;
  double output_sensitivity=0;
  gsl_vector **y,**fs;   // arrays with index structure y(t_j,u_k)=y[c*T+j]
  gsl_matrix **yS,**fyS; //
  gsl_vector *ref_fy;
  gsl_matrix *ref_fyS;

  ref_fy=mp->reference_fy;
  ref_fyS=mp->reference_fyS;
  
  y=mp->y;
  fy=mp->fy;
  yS=mp->yS;
  fyS=mp->fyS;
  for (c=0;c<C;c++){
    for (j=0;j<T;j++) {
      //      gsl_vector_div(y[c*T+j],y[c*T+t0]);
      for (k=0;k<D;k++){
	for (i=0;i<F;i++){
	  output_sensitivity=gsl_matrix_get(fyS[c*T+j],k,i);
	  output_sensitivity-=gsl_matrix_get(fyS[c*T+j],k,i)*gsl_vector_get(fy[c*T+j],k,i)/gsl_vector_get(ref_fy,k,i);
	  output_sensitivity*=gsl_vector_get(mp->exp_x_u,k)/gsl_vector_get(ref_fy,i);
	  gsl_matrix_set(oS[c*T+j],k,i,output_sensitivity);
	}
      }
      gsl_vector_div(fy[c*T+j],ref_fy);
    }
  }
  return GSL_SUCCESS;
}

int normalise_by_timepoint(ode_model_parameters *mp){
  /*  normalisation by one of the points in the time series at t0,
   *   where t0 is the index: t[t0] is the time at which the data is
   *   normalised to 1.
   */
  int c,i,j,k;
  int C=mp->input_u->size1;
  int T=mp->t->size1;
  int F=mp->fy[0]->size2;
  int D=mp->D;
  int N=mp->y[0]->size2;
  double output_sensitivity=0;
  gsl_vector **y,**fs;
  gsl_matrix **yS,**fyS;
  y=mp->y;
  fy=mp->fy;
  yS=mp->yS;
  fyS=mp->fyS;

  t0=gsl_matrix_get(mp->normalisation,0,0);

  for (c=0;c<C;c++){
    for (j=0;j<T;j++){
      for (k=0;k<D;k++){
	for (i=0;i<F;i++){
	  output_sensitivity=gsl_matrix_get(fyS[c*T+j],k,i);
	  output_sensitivity-=gsl_matrix_get(fyS[c*T+j],k,i)*gsl_vector_get(fy[c*T+j],k,i)/gsl_vector_get(fy[c*T+t0],k,i);
	  output_sensitivity*=gsl_vector_get(mp->exp_x_u,k)/gsl_vector_get(fy[c*T+t0],i);
	  gsl_matrix_set(oS[c*T+j],k,i,output_sensitivity);
	}
      }
    }	     
  }
  return GSL_SUCCESS;
}

int normalise_by_state_var(ode_model_parameters *mp){
  /*  normalisation by one of the points in the time series at t0,
   *   where t0 is the index: t[t0] is th etime at which the data is
   *   normalised to 1.
   */
  int c,i,j,k,t0=mp->t0;
  int C=mp->input_u->size1;
  int T=mp->t->size1;
  int F=mp->y[0]->size2;
  int D=mp->D;
  double normalised_fy;
  double output_sensitivity=0;
  gsl_vector **y,**fs;
  gsl_matrix **yS,**fyS;
  gsl_vector_view norm_fy, norm_t;
  gsl_vector *n_fy, *n_t;
  int i_fy, i_t;
    
  norm_fy= gsl_matrix_row(mp->normalisation,0);
  norm_t = gsl_matrix_row(mp->normalisation,1);
  n_fy=&(norm_fy.vector); // output function i is normalised by function n_fy[i]
  n_t =&(norm_t.vector);  // at time index n_t[i]
  
  y=mp->y;
  fy=mp->fy;
  yS=mp->yS;
  fyS=mp->fyS;
  for (c=0;c<C;c++){
    for (j=0;j<T;j++) {
      //      gsl_vector_div(y[c*T+j],y[c*T+t0]);
      for (k=0;k<D;k++){
	for (i=0;i<F;i++){
	  // find the normalising denominator for each sensitivity, by index.
	  i_fy=(int) gsl_vector_get(n_fy,i);
	  i_t =(int) gsl_vector_get(n_t,i);
	  output_sensitivity=gsl_matrix_get(fyS[c*T+j],k,i);
	  output_sensitivity-=gsl_matrix_get(fyS[c*T+j],k,i)*gsl_vector_get(fy[c*T+j],k,i)/gsl_vector_get(fy[c*T+i_t],k,i_fy);
	  output_sensitivity*=gsl_vector_get(mp->exp_x_u,k)/gsl_vector_get(fy[c*T+i_t],i_fy);
	  gsl_matrix_set(oS[c*T+j],k,i,output_sensitivity);
	}
      }
      for (i=0;i<F;i++) {
	// find the normalising denominator for each output function, by index.
	i_fy=(int) gsl_vector_get(n_fy,i);
	i_t =(int) gsl_vector_get(n_t,i);
	normalised_fy=gsl_vector_get(fy[c*T+j],i)/gsl_vector_get(fy[c*T+i_t],i_fy);
	gsl_vector_set(fy[c*T+j],i,normalised_fy);
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
int Posterior(const double* x,  void* model_params, double* fx, 
                     double* dfx, double* FI);


int LikelihoodComplexNorm(ode_model_parameters *mp, double *l, double *dl, double *FI){

  /* Here, the ode integration is done
   * 
   * «l» is a scalar, the return slot of the log-likelihood
   * dl is the parameter gradient of the log likelihood (size=D)
   *
   * cvode returns functions «fy» of the state variables «y».
   * fy=C*y;
   * Note P is the cvode related number of ode parameters
   * while D is the number of MCMC related parameters
   * P=D+U
   * where U is the number of input parameters
   * y[i]: i=0,...,D-1
   * u[i]: i=0,...,U-1
   * cvode sees:
   * p[i]: i=0,...,D-1,D,...,P-1
   *         yyyyyyyyy,uuuuuuuuu
   *
   * FI has to be initialized with zeros(D,D)
   */
  int i,j,k,f,c,u,T,C,F,P,U,N;
  int D=mp->D;
  realtype tout;
  ode_model *model;
  ode_solver *solver;
  double diff;
  gsl_vector *y,*fy;
  gsl_matrix *yS,*fyS;
  gsl_vector *t;
  int iy;

  t=mp->t;
  solver=mp->solver;
  model=solver->odeModel;
  T=t->size;  
  C=mp->input_u->size1;  
  U=mp->input_u->size2;
  P=(mp->D)+U; // P=D+U
  l[0]=0;
  N=ode_model_getN(model);
  F=ode_model_getF(model);
  //  printf("[L] D=%i\tF=%i\tU=%i\tC=%i\tT=%i\tP=%i\n",D,F,U,C,T,P);

  
      
  
  for (c=0; c<C; c++){// loop over different experimental conditions
    // write inputs into the ode parameter vector    
    for (u=0;u<U;u++) gsl_vector_set(mp->exp_x_u,D+u,gsl_matrix_get(mp->input_u,c,u));

    ode_solver_reinit(solver, mp->t0, 0, N,
		      mp->exp_x_u->data,
		      mp->exp_x_u->size);
    
    ode_solver_reinit_sens(solver, mp->yS0->data, P, N);
    for (i=0; i < T; i++){
      y=mp->y[c*T+i];
      fy=mp->fy[c*T+i];
      yS=mp->yS[c*T+i];
      fyS=mp->fyS[c*T+i];

      int CVerror =  ode_solver_solve(solver, gsl_vector_get(t,i), y->data, &tout);
      if (CVerror) {
	fprintf(stderr, "ODE solver failed with ERROR = %i.\n",CVerror);
	// return a rejection message; the Likelihood is not defined for this argument;
	return GSL_EDOM;
      }
      // get sensitivities and output function values for the calculated ODE solutions
      if (ode_model_has_funcs(model)) {
	ode_solver_get_func(solver, tout, y->data, fy->data);
      }
      else {
	error("ode model has no functions");
      }
      
      if (ode_model_has_sens(model)){
	ode_solver_get_sens(solver, tout, yS->data);
	if (ode_model_has_funcs_sens(model)){
	  ode_solver_get_func_sens(solver, tout, y->data, yS->data, fyS->data);
	} // end if has func sens
      } // end if has sens
      else{
	error("ode model has no sensitivities.");
      }	
      /* normalise ode solution to reflect data
       * normalisation. Normalisation will be performed "in place",
       * i.e. the ODE solution will be overwritten
       */
      switch (mp->normalisation_method){
      case DATA_NORMALISED_BY_TIMEPOINT: break;
      case DATA_NORMALISED_BY_REFERENCE: break;
      case DATA_NORMALISED_BY_STATE_VAR: break;
      default:
      }
    } // end for loop for time points
  } //end for different experimental conditions (i.e. inputs)
  return GSL_SUCCESS;
}

int Likelihood(ode_model_parameters *mp, double *l, double *dl, double *FI){

  /* Here, the ode integration is done
   * 
   * «l» is a scalar, the return slot of the log-likelihood
   * dl is the parameter gradient of the log likelihood (size=D)
   *
   * cvode returns functions «fy» of the state variables «y».
   * fy=C*y;
   * Note P is the cvode related number of ode parameters
   * while D is the number of MCMC related parameters
   * P=D+U
   * where U is the number of input parameters
   * y[i]: i=0,...,D-1
   * u[i]: i=0,...,U-1
   * cvode sees:
   * p[i]: i=0,...,D-1,D,...,P-1
   *         yyyyyyyyy,uuuuuuuuu
   *
   * FI has to be initialized with zeros(D,D)
   */
  int i,j,k,f,c,u,T,C,F,P,U,N;
  int D=mp->D;
  realtype tout;
  ode_model *model;
  ode_solver *solver;
  double diff;
  gsl_vector *y,*fy,*yS,*fyS;
  gsl_vector_view ref_fy_ti;
  gsl_vector *rfyi,*t;
  int iy;

  t=mp->t;
  solver=mp->solver;
  model=solver->odeModel;
  T=t->size;  
  C=mp->input_u->size1;  
  U=mp->input_u->size2;
  P=(mp->D)+U; // P=D+U
  l[0]=0;
  N=ode_model_getN(model);
  F=ode_model_getF(model);
  //  printf("[L] D=%i\tF=%i\tU=%i\tC=%i\tT=%i\tP=%i\n",D,F,U,C,T,P);

  for (c=-1; c<C; c++){// loop over different experimental conditions
    // write inputs into the ode parameter vector    
    if (c<0){ // reference experiment
      for (u=0;u<U;u++) gsl_vector_set(mp->exp_x_u,D+u,gsl_vector_get(mp->reference_u,u));
    }else{ 
      for (u=0;u<U;u++) gsl_vector_set(mp->exp_x_u,D+u,gsl_matrix_get(mp->input_u,c,u));
    }

    ode_solver_reinit(solver, mp->t0, 0, N,
		      mp->exp_x_u->data,
		      mp->exp_x_u->size);
    
    ode_solver_reinit_sens(solver, mp->yS0->data, P, N);

    for (i=0; i < t->size; i++){
      if (c<0){ // reference experiment
	y=&(mp->reference_y->data[i*N]);        // access the correct block for
	fy=&(mp->reference_fy->data[i*F]);      // the current time point t[i]
	yS=&(mp->reference_yS->data[i*P*N]);    // so for fyS_{ik}^j the general
	fyS=&(mp->reference_fyS->data[i*P*F]);  // index is: i*P*F+j*F+k
      }else{
	y=mp->y->data;
	fy=mp->fy->data;
	yS=mp->yS->data;
	fyS=mp->fyS->data;
      }

      if ((c<0) && !(mp->data_is_relative) ){
	/* set reference_fyS to 0.0 and reference_fy to 1.0, i.e. the
	   reference is a constant 1.0 */
	gsl_matrix_set_zero(mp->reference_fyS);
	gsl_matrix_set_all(mp->reference_fy,1.0);
      }else{
	int CVerror =  ode_solver_solve(solver, gsl_vector_get(t,i), y, &tout);
	//    printf("[Likelihood (i=%i, c=%i) t=%g] y = ",i,c,gsl_vector_get(t,i));
        //for (iy=0;iy<N;iy++) printf(" %g ",y[iy]); printf("\n");

	if (CVerror) {
	  fprintf(stderr, "ODE solver failed with ERROR = %i.\n",CVerror);
	  // return a rejection message; the Likelihood is not defined for this argument;
	  return GSL_EDOM;
	}
	//printf("ode_model_has_functions: %i\n",ode_model_has_funcs(model));
	if (ode_model_has_funcs(model)) {
	  ode_solver_get_func(solver, tout, y, fy);
	}
	else {
	  error("ode model has no functions");
	}

	if (ode_model_has_sens(model)){
	  ode_solver_get_sens(solver, tout, yS);
	  if (ode_model_has_funcs_sens(model)){
	    ode_solver_get_func_sens(solver, tout, y, yS, fyS);
	  } // end if has func sens
	} // end if has sens
	else{
	  error("ode model has no sensitivities.");
	}	
      }	
      if (c>=0){ /* only start calculating the likelihood, gradient
		  * and FI when reference observations are done
		  */
	ref_fy_ti=gsl_matrix_row(mp->reference_fy,i); // point to the right row
	rfyi=&(ref_fy_ti.vector);
	gsl_vector_div(mp->fy,rfyi);
	for (f=0;f<F;f++) {
	  diff=gsl_matrix_get(mp->Data,c*T+i,f)-gsl_vector_get(mp->fy,f);
	  diff/=gsl_matrix_get(mp->sdData,c*T+i,f);
	  l[0]+=-0.5*gsl_pow_2(diff);
	} // end likelihood sum
	/* printf("l = %g\n",l[0]); */
	/* gsl_printf("exp_x_u",mp->exp_x_u,0); */
	/* gsl_printf("fyS",mp->fyS,1); */
	/* gsl_printf("fy",mp->fy,0); */
	/* gsl_printf("reference_fyS",mp->reference_fyS,1); */
	for (j=0;j<D;j++) {
	  /*
	   * the dense index of sensitivities is the state (or
	   * "function") index 
	   * yS[i*N+j] = dy[j]/dp[i];
	   * fyS[i*F+j] = dfy[j]/dp[i];
	   */
	  for (f=0;f<F;f++){ // calculate observational sensitivities, i.e. ~{fy/ref_fy}S 
	    gsl_matrix_set(mp->oS,j,f,
			   gsl_vector_get(mp->exp_x_u,j)*
			   (gsl_matrix_get(mp->fyS,j,f)
			    -gsl_vector_get(mp->fy,f)*gsl_matrix_get(mp->reference_fyS,i*P+j,f))
			   /gsl_vector_get(rfyi,f)
			   );
	  }
	  if (dl!=NULL){
	    dl[j]=0;
	    for (f=0;f<F;f++){ // loop over y functions
	      diff=gsl_matrix_get(mp->Data,c*T+i,f)-gsl_vector_get(mp->fy,f);
	      diff/=gsl_pow_2(gsl_matrix_get(mp->sdData,c*T+i,f));
	      dl[j]+=diff*gsl_matrix_get(mp->oS,j,f);
	    } // end loop over functions
	    //printf("dl[%i] = %g \n",j,dl[j]);
	  } // end if dl != NULL
	} // end gradient loop
	
	/* calculate Fisher Information; we are still inside the time loop
	   (over t[i])*/
	if (FI!=NULL){	  
	  for (j=0;j<D;j++){ // loop over rows of FI, 0,...,D-1
	    for (k=0;k<D;k++){ // loop over columns of FI, 0,...,D-1
	      for (f=0;f<F;f++){ // loop over functions
		FI[j*D+k]+=gsl_matrix_get(mp->oS,j,f)*gsl_matrix_get(mp->oS,k,f)
		  /gsl_pow_2(gsl_matrix_get(mp->sdData,c*T+i,f)); // squared
	      }// end loop over functions
	    }// end loop over columns
	  }// end loop over rows	  
	} // end if FI!=NULL	
      } // end if c>=0
    } // end for loop for time points
  } //end for different experimental conditions (i.e. inputs)
  return GSL_SUCCESS;
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
  printf("\t\t\tresume from last sampled MCMC point. Only the last MCMC position is read from the file named «resume.double». Everything else can be changed.")
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


int main (int argc, char* argv[]) {
  int D = 0; // number of MCMC sampling variables, i.e. model parameters
  int C = 0; // number of experimental conditions, i.e different input vectors
  int i,j;
  char *cfilename=NULL;
  char lib_name[128];
  ode_model_parameters omp;
  FILE *cnf;   // configuration file, with file name: cfilename
  char sample_file[128]="sample.dat"; // filename for sample output
  char *x_sample_file=NULL; // filename for sample output x(t,p)
  char *y_sample_file=NULL; // filename for sample output y(t,p)
  FILE *oFile; // will be the file named «sample_file»
  FILE *rFile; // last sampled value will be written to this file
  char resume_filename[128]="resume.double";
  int output_is_binary=0;
  double seed = 1;
  int Tuning = 1;
  int sampling_action=SMPL_FRESH;
  size_t resume_count;
  gsl_error_handler_t *gsl_error_handler;
  
  gsl_error_handler = gsl_set_error_handler_off();

  //  strcpy(sample_file,"sample.dat");
  main_options cnf_options;
  cnf_options.output_file=sample_file;
  cnf_options.library_file=lib_name;
  cnf_options.target_acceptance=-0.5;
  cnf_options.initial_stepsize =-0.1; // not used here in smmala
  cnf_options.sample_size=-10;

  for (i=0;i<argc;i++){
    if (strcmp(argv[i],"-c")==0) cfilename=argv[i+1];
    else if (strcmp(argv[i],"--resume")==0 || strcmp(argv[i],"-r")==0) sampling_action=SMPL_RESUME;
    else if (strcmp(argv[i],"-l")==0) strcpy(cnf_options.library_file,argv[i+1]);
    else if (strcmp(argv[i],"-n")==0) Tuning=0;
    else if (strcmp(argv[i],"-s")==0) cnf_options.sample_size=strtol(argv[i+1],NULL,0);
    else if (strcmp(argv[i],"-o")==0) strcpy(cnf_options.output_file,argv[i+1]);
    else if (strcmp(argv[i],"-b")==0) output_is_binary=1;
    else if (strcmp(argv[i],"-a")==0) cnf_options.target_acceptance=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-i")==0) cnf_options.initial_stepsize=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"--seed")==0) seed=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-h")==0) {print_help(); exit(0);}
    
    //else printf("unknown option {%s}.\n",argv[i]);
  }

  /* load model */
  ode_model *odeModel = ode_model_loadFromFile(lib_name);  /* alloc */
  if (odeModel == 0) {
    fprintf(stderr, "# Library %s could not be loaded.\n",lib_name);
    exit(1);
  } else printf( "# Library %s loaded.\n",lib_name);

  ode_solver* solver = ode_solver_alloc(odeModel); /* alloc */
  if (solver == 0) {
    fprintf(stderr, "# Solver %s could not be created.\n",lib_name);
    ode_model_free(odeModel);
    exit(1);
  } else printf("# Solver for %s created.\n",lib_name);

  /* init solver */
  realtype solver_param[3] = {ODE_SOLVER_ABS_ERR, ODE_SOLVER_REL_ERR, 0};
  
  /* define local variables for parameters and inital conditions */
  int N = ode_model_getN(odeModel);
  int P = ode_model_getP(odeModel);
  int F = ode_model_getF(odeModel);
  problem_size ps;
  ps.N=N;
  ps.P=P;
  ps.F=F;
  double y[N];
  ode_model_get_initial_conditions(odeModel, y, N);
	
  double p[P];
  ode_model_get_default_params(odeModel, p, P);
	
  //  double t_old = solver_param[2];

  omp.solver=solver;

  cnf=fopen(cfilename,"r");
  if (cnf!=NULL){
    printf("# reading configuration.\n");
    parse_config(cnf,&omp,&ps,&cnf_options);
  } else {
    fprintf(stderr,"# Could not open config file %s.\n",cfilename);
    exit(1);
  }
  fclose(cnf);
 
  cnf_options.initial_stepsize=fabs(cnf_options.initial_stepsize);
  cnf_options.target_acceptance=fabs(cnf_options.target_acceptance);
  cnf_options.sample_size=fabs(cnf_options.sample_size);

  /* gsl_printf("data",omp.Data,1); */
  /* gsl_printf("sd data",omp.sdData,1); */
  /* gsl_printf("u",omp.input_u,1); */
  /* gsl_printf("t",omp.t,0); */
  
  printf("# init ivp: t0=%g\n",omp.t0);
  ode_solver_init(solver, omp.t0, y, N, p, P);
  ode_solver_setErrTol(solver, solver_param[1], &solver_param[0], 1);
  if (ode_model_has_sens(odeModel)) {
    ode_solver_init_sens(solver, omp.yS0->data, P, N);
    printf("# sensitivity analysis initiated.\n");
  }

  printf("# allocating memory for the SMMALA model.\n");
	
  smmala_model* model = smmala_model_alloc(Posterior, NULL, &omp); // ode_model_parameters
  /* initial parameter values */
  double init_x[omp.D];
	
  /* allocate a new RMHMC MCMC kernel */
  printf("# allocating memory for a new SMMALA MCMC kernel.\n");
  /*mcmc_kernel* smmala_kernel_alloc(
    int N, 
    double step_size,
    smmala_model* model_function, 
    unsigned long int seed)
  */
  mcmc_kernel* kernel = smmala_kernel_alloc(omp.D, cnf_options.initial_stepsize, model, seed, cnf_options.target_acceptance);
  printf("# initializing MCMC.\n");

  /* initialise MCMC */
  for (i=0;i<omp.D;i++) init_x[i]=gsl_vector_get(omp.prior_mu,i);
  if (sampling_action==SMPL_RESUME){
    rFile=fopen(resume_filename,"r");
    if (rFile==NULL) {
      printf("Could not open resume file. Starting from: ");
      gsl_printf("prior mean" omp.prior_mu,0);
    } else {
      size_t resume_count=fread(init_x, sizof(double), D, rFile);
      fclose(rFile);
      if (resume_count!=D) error("Reading from resume file returned a wrong number of values.");
    }
  }
  mcmc_init(kernel, init_x);

  printf("# test evaluation of Posterior function.\n");
  /*inits*/
  double fx;
  double dfx[omp.D];
  double FI[omp.D*omp.D];
  Posterior(init_x, &omp, &fx, dfx, FI);
  printf("# x=ones(%i); Posterior(x)=%f;\n",omp.D,fx);

  /* print first sample, initial values in init_x */
  mcmc_print_sample(kernel, stdout);

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


  size_t BurnInSamples = 7 * (int) sqrt(cnf_options.sample_size);
  if (sampling_action==SMPL_FRESH){
    fprintf(stdout, "# Burn-in phase.\n");
    /* "find mode" Loop */
    for (it = 0; it < BurnInSamples; it++) {
      omp.beta=((double) it)/((double) BurnInSamples);
      omp.beta=gsl_pow_6(omp.beta);
      omp.beta=omp.beta/(gsl_pow_6(0.5)+omp.beta);
      mcmc_sample(kernel, &acc);
      acc_c += acc;
      /* print sample */
      mcmc_print_sample(kernel, stdout);
		
      /* Addapt MCMC parameters every 100 burn-in samples */
      if ( ((it + 1) % 100) == 0 ) {
	acc_rate = (double)acc_c / (double)100;
	fprintf(stdout, "# Iteration: %li\tAcceptance rate: %.2g\t",it, acc_rate);
	/* print log and statistics */
	mcmc_print_stats(kernel, stdout);
	/* adapt parameters */
	mcmc_adapt(kernel, acc_rate);
	acc_c = 0;
      }
    }

    omp.beta=1.0;

    /* Burn In Loop and Tuning*/
    for (it = 0; it < BurnInSamples; it++) {
      mcmc_sample(kernel, &acc);
      acc_c += acc;
      /* print sample */
      mcmc_print_sample(kernel, stdout);
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
  fclose(rFile)
  ct=clock()-ct;
  fclose(oFile);
  printf("# computation time spend sampling: %f s\n",((double) ct)/((double) CLOCKS_PER_SEC));


  /* clear memory */
  smmala_model_free(model);
  mcmc_free(kernel);
  ode_model_parameters_free(&omp);
  return EXIT_SUCCESS;
}



/* Updates the mvNormParams structure with the inverse of the covariance matrix.
 */
/*void ParamsInvertCovariance(ode_model_parameters* params){
//	 do it only once since the variance does not depent on x 
	int D = params->D;
	if (params->init == 0) {
		gsl_matrix_view invV_v = gsl_matrix_view_array(params->Precision, D, D);
		gsl_matrix_const_view V_v = gsl_matrix_const_view_array(params->Variance, D, D);
		
		gsl_matrix_memcpy(&invV_v.matrix, &V_v.matrix);
		gsl_linalg_cholesky_decomp(&invV_v.matrix);
		gsl_linalg_cholesky_invert(&invV_v.matrix);
		
		params->init = 1;
	}
}
*/

/* Calculates the unormalised log? posterior fx, the gradient dfx, the Fisher information FI 
 * and the partial derivatives of the Fisher information matrix for ode models
 */
int Posterior(const double* x,  void* model_params, double* fx, double* dfx, double* FI){
	
  ode_model_parameters* params =  (ode_model_parameters*) model_params;
  double prior_value=0;            // the prior distribution related sum of squares
  int i,j,F,T,D,U;
  gsl_matrix *Data,*sd;
  gsl_vector *fy;
  gsl_matrix *inv_cov;
  gsl_vector *prior_diff, *prior_mv;
  double l;
  int logL_stat;

  U=params->input_u->size2;
  D=params->D;
  prior_diff=params->prior_tmp_a;
  prior_mv=params->prior_tmp_b;

  fy=params->fy;
  Data=params->Data;
  sd=params->sdData;

  inv_cov=params->prior_inverse_cov;

  T=(params->t)->size;
  F=(params->fy)->size;
  // init FI and dFI
  fx[0]=0;
  for (i=0;i<D;i++) gsl_vector_set(params->exp_x_u,i,gsl_sf_exp(x[i]));
  for (i=0;i<D*D;i++) FI[i]=0;
 
  logL_stat=Likelihood(params, fx, dfx, FI);
  //  printf("Likelihood: fx = %f\n",*fx);
  //  printf("Likelihood: dfx = ["); for (i=0;i<D;i++) printf("%f\t",dfx[i]); printf("]\n");
  //  printf("[l FI]\n");
  /*for (i=0;i<D;i++){
    for (j=0;j<D;j++) printf("%f\t",FI[i*D+j]); printf("\n");
    }*/
  //printf("[/l FI]\n");

  for (i=0;i<D;i++) gsl_vector_set(prior_diff,i,x[i]-gsl_vector_get(params->prior_mu,i));

  gsl_blas_dsymv(CblasUpper,-0.5,inv_cov,prior_diff,0.0,prior_mv);
  gsl_blas_ddot(prior_diff,prior_mv,&prior_value);
  //  printf("PosteriorFI: prior_value=%f\n",prior_value);
  fx[0]*=params->beta;
  fx[0]+=prior_value;
  /* here the "sum_of_squares" belongs to the log-likelihood
   * while the prior_value is the prior-related part of the log-posterior
   * these values omit all constants.
   */

  // now for the gradient:
  for (i=0; i < D ; i++){
    dfx[i]*= params->beta;
    dfx[i]+= gsl_vector_get(prior_mv,i); // prior component
  }
  
  /* FI is constant, it does not depend on the parameters x ... */
  for (i=0; i < D*D ; i++){
    FI[i] *= gsl_pow_2(params->beta);
    FI[i] += params->prior_inverse_cov->data[i];
  }
  /* printf("[p FI]\n"); */
  /* for (i=0;i<D;i++){ */
  /*   for (j=0;j<D;j++) printf("%f\t",FI[i*D+j]); printf("\n"); */
  /* } */
  /* printf("[/p FI]\n"); */

  return logL_stat;
}

