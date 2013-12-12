/*
 *  ode_rmhmc.c
 *  
 *  based on multivariate normal implementation: mvNorm_rmhmc.c
 *  modified by Andrei Kramer 
 *  andrei.kramer@ist.uni-stuttgart.de
 *
 *  Copyright (C) 2012 Vassilios Stathopoulos stathv@gmail.com
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
#include <gsl/gsl_errno.h>
#include "read_cnf.h"
#include "../mcmc/RMHMC.h"
#include "../ode/ode_model.h"

/* Auxiliary structure with working storage and aditional parameters for
 * a multivariate normal model with known covariance matrix and zero mean.
 *
 */



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
                     double* dfx, double* FI, double* dFI );


/* Calculates the Fisher information matrix for a multivariate normal.
 * Arguments:
 *  x:				values for the model parameters.
 *	model_params:	a pointer to mvNormParams struct.
 * Returns:
 *  FI:				the Fisher information evaluated at x. A double array with size N*N whith 
 *					elements stored in row major order. The i,j element of the FI is stored at FI[i*N + j].
 * The function should return 0 for success or a non-zero value indicating failure.
 */
int PosteriorFI(const double* x,  void* model_params, double* FI);


/* Auxiliary function.
 * Updates the mvNormParams structure with the inverse of the covariance matrix.
 */
void ParamsInvertCovariance(ode_model_parameters* params);

int gsl_printf(char *name, void *gsl_thing, int is_matrix){
  gsl_matrix *A;
  gsl_vector *x;
  
  printf("[%s]\n",name);
  if (is_matrix){
    A=(gsl_matrix*) gsl_thing;
    gsl_matrix_fprintf (stdout, A, "%g");
  }else{
    x=(gsl_vector*) gsl_thing;
    gsl_vector_fprintf (stdout, x, "%g");
  }
  printf("[/%s]\n",name);
  return GSL_SUCCESS;
}

int Likelihood(ode_model_parameters *mp, double *l, double *dl, double *FI, double *dFI){

  /* Here, the ode integration is done
   * 
   * «l» is a scalar, the return slot of the log-likelihood
   * dl is the parameter gradient of the log likelihood (size=P)
   *
   * cvode returns functions «fy» of the state variables «y».
   * 
   * Note P is the cvode related number of ode parameters
   * while D is the number of MCMC related parameters
   * P=D+U
   * where U is the number of input parameters
   * x[i]: i=0,...,D-1
   * u[i]: i=0,...,U-1
   * cvode sees:
   * p[i]: i=0,...,D-1,D,...,P-1
   *         xxxxxxxxx,uuuuuuuuu
   *
   * FI has to be initialized with zeros(D,D)
   */
  int it,i,j,k,f,c,u,T,C,F,P,U,N;
  int cvodes_N,cvodes_F;
  int D=mp->D;
  realtype tout;
  ode_model *model;
  ode_solver *solver;
  double diff;
  double *y,*yS,*ry,*ryS;
  double *dyS,*dryS;
  // to be able to use all kinds of gsl functions we create some views
  gsl_vector_view y_vv,ry_vv;
  gsl_matrix_view yS_mv,ryS_mv;
  gsl_matrix_view dyS_mv,dryS_mv;
  gsl_vector *t;

  t=mp->t;
  solver=mp->solver;
  model=solver->odeModel;
  T=t->size;  
  C=mp->input_u->size1;  
  U=mp->input_u->size2;
  P=(mp->D)+U; // P=D+U
  l[0]=0;
  N=mp->output_C->size2;
  F=mp->output_C->size1;

  cvodes_N=ode_model_getN(model);
  cvodes_F=ode_model_getF(model);
  //  printf("cvodes_N=%i\n",cvodes_N);
  //  printf("[L] D=%i\tF=%i\tU=%i\tC=%i\tT=%i\tP=%i\n",D,F,U,C,T,P);

  for (c=-1; c<C; c++){// loop over different experimental conditions
    // write inputs into the ode parameter vector    
    if (c<0){ // reference experiment
      for (u=0;u<U;u++) gsl_vector_set(mp->exp_x_u,D+u,gsl_vector_get(mp->reference_u,u));
    }else{ 
      for (u=0;u<U;u++) gsl_vector_set(mp->exp_x_u,D+u,gsl_matrix_get(mp->input_u,c,u));
    }
    ode_solver_reinit(solver,	mp->t0,	0, cvodes_N,
		      mp->exp_x_u->data,
		      mp->exp_x_u->size);
    
    for (it=0; it < t->size; it++){
      if (c<0){ // reference experiment
	y=&(mp->cvodes_reference_state[it*cvodes_N]);
      }else{
	y=mp->cvodes_state;
      }

      if ((c<0) && !(mp->data_is_relative) ){
	/* set reference_fyS to 0.0 and reference_fy to 1.0, i.e. the
	   reference is a constant 1.0 */
	gsl_matrix_set_zero(mp->drfyS);
	gsl_matrix_set_zero(mp->rfyS);
	gsl_vector_set_all(mp->rfy,1.0);
      }else{
	int CVerror =  ode_solver_solve(solver, gsl_vector_get(t,it), y, &tout);
	if (CVerror) {
	  fprintf(stderr, "ODE solver failed.\n");
	  exit(2);
	}
      }	
      if (c>=0){ /* only start calculating the likelihood, gradient
		  * and FI when reference observations are done
		  */
	yS=&(y[N]);      // yS comes after the state variables
	dyS=&(yS[P*N]);  // dyS comes after the sensitivity

	// set up some views to use gsl
	y_vv=gsl_vector_view_array(y,N); 
	yS_mv=gsl_matrix_view_array(yS,P,N);
	dyS_mv=gsl_matrix_view_array(dyS,P*P,N);

	// fy=C*y; fyS=C*yS
	gsl_blas_dgemv(CblasNoTrans,1.0,mp->output_C,&(y_vv.vector),0.0,mp->fy);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &(yS_mv.matrix), mp->output_C, 0.0, mp->fyS);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &(dyS_mv.matrix), mp->output_C, 0.0, mp->dfyS);

	if (mp->data_is_relative){
	  ry=&(mp->cvodes_reference_state[it*cvodes_N]);
	  ryS=&(ry[N]);
	  dryS=&(ryS[P*N]);

	  ry_vv=gsl_vector_view_array(ry,N); 
	  ryS_mv=gsl_matrix_view_array(ryS,P,N);
	  dryS_mv=gsl_matrix_view_array(dryS,P*P,N);
	  
	  // rfy=C*ry; rfyS=C*ryS
	  gsl_blas_dgemv(CblasNoTrans, 1.0,mp->output_C,&(ry_vv.vector), 0.0,mp->rfy);
	  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &(ryS_mv.matrix), mp->output_C, 0.0, mp->rfyS);
	  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &(dryS_mv.matrix), mp->output_C, 0.0, mp->drfyS);
	  // fy:=fy/rfy;
	  gsl_vector_div(mp->fy,mp->rfy);
	}
	
	// from this point on fy is a relative value if applicable:
	for (f=0;f<F;f++) {
	  diff=gsl_matrix_get(mp->Data,c*T+it,f)-gsl_vector_get(mp->fy,f);
	  diff/=gsl_matrix_get(mp->sdData,c*T+it,f);
	  l[0]+=-0.5*gsl_pow_2(diff);
	} // end likelihood sum

	for (j=0;j<D;j++) {
	  /*
	   * the dense index of sensitivities is the state (or
	   * "function") index 
	   * yS[i*N+j] = dy[j]/dp[i];
	   * fyS[i*F+j] = dfy[j]/dp[i];
	   *
	   * the index structur of sensitivity gradients: d{any}S
	   * dfyS[j*P*F+k*F+f]=d fyS[j*F+f]/dp[k]; 
	   *
	   * => k is more dense (faster) than j, which is
	   * different in dFI[k*D*D+i*D+j];
	   */
	  for (f=0;f<F;f++){ // calculate observational sensitivities, i.e. ~{fy/ref_fy}S 
	    gsl_matrix_set(mp->oS,j,f,
			   (gsl_matrix_get(mp->fyS,j,f)
			    -gsl_vector_get(mp->fy,f)*gsl_matrix_get(mp->rfyS,j,f))
			   /gsl_vector_get(mp->rfy,f)
			   );
	  }
	  for (k=0;k<D;k++){
	    for (f=0;f<F;f++){ 
	      // calculate the derivatives of the
	      // observational sensitivities,
	      // i.e. ~d{fy/ref_fy}S
	      gsl_matrix_set(mp->doS,j*D+k,f,
			     (gsl_matrix_get(mp->dfyS,j*P+k,f) 
			      - gsl_matrix_get(mp->oS,k,f)*gsl_matrix_get(mp->rfyS,j,f)
			      - gsl_vector_get(mp->fy,f)*gsl_matrix_get(mp->drfyS,j*P+k,f))
			     /gsl_vector_get(mp->rfy,f)
			     + (gsl_vector_get(mp->fy,f)*gsl_matrix_get(mp->rfyS,j,f)
				-gsl_matrix_get(mp->fyS,j,f))*gsl_matrix_get(mp->rfyS,k,f)
			     /gsl_pow_2(gsl_vector_get(mp->rfy,f))
			     );
	      mp->doS->data[j*D*F+k*F+f]*=gsl_vector_get(mp->exp_x_u,j)*gsl_vector_get(mp->exp_x_u,k);
	    }//end f loop
	  }//end k loop
	  // this involves a δ_jk, therefore this expression is outside of the k loop, since k=j
	  for (f=0;f<F;f++){ 
	    mp->doS->data[j*D*F+j*F+f]+=gsl_matrix_get(mp->oS,j,f)*gsl_vector_get(mp->exp_x_u,j);
	    mp->oS->data[j*F+f]*=gsl_vector_get(mp->exp_x_u,j);
	  }

	  if (dl!=NULL){
	    dl[j]=0;
	    for (f=0;f<F;f++){ // loop over y functions
	      diff=gsl_matrix_get(mp->Data,c*T+it,f)-gsl_vector_get(mp->fy,f);
	      diff/=gsl_pow_2(gsl_matrix_get(mp->sdData,c*T+it,f));
	      dl[j]+=diff*gsl_vector_get(mp->exp_x_u,j)
		*gsl_matrix_get(mp->oS,j,f);
	    } // end loop over functions
	  } // end if dl != NULL
	} // end gradient loop
	
	/* calculate Fisher Information; we are still inside the time loop
	   (over t[i])*/
	if (FI!=NULL){
	  for (j=0;j<D;j++){ // loop over rows of FI, 0,...,D-1
	    for (k=0;k<D;k++){ // loop over columns of FI, 0,...,D-1
	      for (f=0;f<F;f++){ // loop over functions
		FI[j*D+k]+=gsl_matrix_get(mp->oS,j,f)*gsl_matrix_get(mp->oS,k,f)
		  /gsl_pow_2(gsl_matrix_get(mp->sdData,c*T+it,f));
	      }// end loop over functions
	    }// end loop over columns
	  }// end loop over rows
	} // end if FI!=NULL
	if (dFI!=NULL){
	  for (k=0;k<D;k++){
	    for (i=0;i<D;i++){
	      for (j=0;j<D;j++){
		for (f=0;f<F;f++){
		  dFI[k*D*D+i*D+j]+=(gsl_matrix_get(mp->doS,i*D+k,f)*gsl_matrix_get(mp->oS,j,f)
				     +gsl_matrix_get(mp->oS,i,f)*gsl_matrix_get(mp->doS,j*D+k,f))
		                      /gsl_pow_2(gsl_matrix_get(mp->sdData,c*T+it,f));
		}
	      }
	    }
	  }
	} // end if dFI!=NULL	
      } // end if c>=0
    } // end for loop for time points
  } //end for different experimental conditions (i.e. inputs)

  //print all calculated values:
  /*
  gsl_printf("C",mp->output_C,1);
  gsl_printf("fy",mp->fy,0);
  gsl_printf("fyS",mp->fyS,1);
  gsl_printf("dfyS",mp->dfyS,1);
  */
  //lets print the derivative of the sensitivity:
  /*
  printf("e^x[0]=%g\n",gsl_vector_get(mp->exp_x_u,0));
  gsl_printf("oS",mp->oS,1);
  gsl_printf("doS",mp->doS,1); 
  */
  //print FI
  /*
  printf("[FI]\n");
    for (i=0;i<D;i++){
      for (j=0;j<D;j++) printf("%f\t",FI[i*D+j]);
      printf("\n");
    }
  printf("[/FI]\n");
  */

  //print dFI
  /*
  printf("[dFI]\n");
  for (k=0;k<D;k++){
    for (i=0;i<D;i++){
      for (j=0;j<D;j++) printf("%f\t",dFI[k*D*D+i*D+j]);
      printf("\n");
    }
    printf("\n");
  }
  printf("[/dFI]\n");
  exit(0);
  */
  return EXIT_SUCCESS;
}

int print_help(){
  printf("Usage:\n");
  printf("-c ./data.cfg\n");
  printf("\t\t\tdata.cfg contains the data points and the conditions of measurement.\n\n");
  printf("-l ./ode_model.so\n");
  printf("\t\t\tode_model.so is a shared library containing the CVODE functions of the model.\n\n");
  printf("-s $N\n");
  printf("\t\t\t$N sample size. default N=10.\n\n");
  printf("-o ./output_file\n");
  printf("\t\t\tfile for output. Output can be binary.\n\n");
  printf("-b\n");
  printf("\t\t\toutput mode: binary.\n\n");
  printf("-i ${step_size}\n");
  printf("\t\t\tsets the initial step size.\n\n");
  printf("--seed\n");
  printf("\t\t\tsets the seed for the random number generator.\n\n");
  return EXIT_SUCCESS;
}



int main (int argc, char* argv[]) {
  // later this should be read from a config file	
  int i,j,k;
  char *cfilename=NULL;
  char lib_name[256];
  ode_model_parameters omp;
  FILE *cnf; // configuration file, with file name: cfilename
  char sample_file[256]="ode_parameter_sample.dat";
  FILE *oFile; // will be the file named «sample_file»
  int output_is_binary=0;

  double leap_frog_steps = 10;
  double fixed_point = 3;
  double seed = 1;
  
  main_options cnf_options;
  cnf_options.output_file=sample_file;
  cnf_options.library_file=lib_name;
  cnf_options.target_acceptance=-0.95; // default values are set negative
  cnf_options.initial_stepsize =-0.01; // this makes it easy to detect, if they were set on the command line 
  cnf_options.sample_size=-10;         // (then: value > 0)

  gsl_set_error_handler_off();

  for (i=0;i<argc;i++){
    if (strcmp(argv[i],"-c")==0) cfilename=argv[i+1];
    else if (strcmp(argv[i],"-l")==0) strcpy(cnf_options.library_file,argv[i+1]);
    else if (strcmp(argv[i],"-s")==0) cnf_options.sample_size=strtol(argv[i+1],NULL,0);
    else if (strcmp(argv[i],"-o")==0) strcpy(cnf_options.output_file,argv[i+1]);
    else if (strcmp(argv[i],"-b")==0) output_is_binary=1;
    else if (strcmp(argv[i],"-i")==0) cnf_options.initial_stepsize=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-a")==0) cnf_options.target_acceptance=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"--seed")==0) seed=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-h")==0) {print_help(); exit(0);}
  }



  /* load model */
  ode_model *odeModel = ode_model_loadFromFile(lib_name);  /* alloc */
  if (odeModel == 0) {
    fprintf(stderr, "Library %s could not be loaded.\n",lib_name);
    exit(1);
  } else printf( "Library %s loaded.\n",lib_name);

  ode_solver* solver = ode_solver_alloc(odeModel); /* alloc */
  if (solver == 0) {
    fprintf(stderr, "Solver %s could not be created.\n",lib_name);
    ode_model_free(odeModel);
    exit(1);
  } else printf("Solver for %s created.\n",lib_name);

  /* init solver */
  realtype solver_param[3] = {ODE_SOLVER_ABS_ERR, ODE_SOLVER_REL_ERR, 0};
  
  /* define local variables for parameters and inital conditions */
  int cvodes_N = ode_model_getN(odeModel);
  int P = ode_model_getP(odeModel);
  int F = ode_model_getF(odeModel);
	
  double y[cvodes_N];
  ode_model_get_initial_conditions(odeModel, y, cvodes_N);
	
  double p[P];
  ode_model_get_default_params(odeModel, p, P);
	
  double t_old = solver_param[2];

  omp.solver=solver;
  // read config
  cnf=fopen(cfilename,"r");
  if (cnf!=NULL){
    printf("reading configuration.\n");
    parse_config(cnf,&omp,&cnf_options);
  } else {
    fprintf(stderr,"Could not open config file %s.\n",cfilename);
    exit(1);
  }
  fclose(cnf);
  
  cnf_options.initial_stepsize=fabs(cnf_options.initial_stepsize);
  cnf_options.target_acceptance=fabs(cnf_options.target_acceptance);
  cnf_options.sample_size=fabs(cnf_options.sample_size);


  printf("Data (%i×%i)\n",omp.Data->size1,omp.Data->size2);

  gsl_printf("data",omp.Data,1);
  gsl_printf("sd data",omp.sdData,1);


  ode_solver_init(solver, t_old, y, cvodes_N, p, P);
  ode_solver_setErrTol(solver, solver_param[1], &solver_param[0], 1);


  printf("allocating memory for the RMHMC model.\n");

  /*  allocate new RMHMC model structure with the user defined functions */
	
  rmhmc_model* model = rmhmc_model_alloc(Posterior, 
					 PosteriorFI, 
					 NULL, &omp); // ode_model_parameters
  /* initial parameter values */
  double init_x[omp.D];
	
  /* allocate a new RMHMC MCMC kernel */
  printf("allocating memory for a new RMHMC MCMC kernel.\n");
  mcmc_kernel* kernel = rmhmc_kernel_alloc(omp.D, cnf_options.initial_stepsize, leap_frog_steps, 
					   fixed_point, model, seed, cnf_options.target_acceptance);
  printf("initialize MCMC.\n");
  /* initialise MCMC */
  for (i=0;i<omp.D;i++) init_x[i]=gsl_vector_get(omp.prior_mu,i);
  mcmc_init(kernel, init_x);

  printf("test evaluation of Posterior function.\n");
  /*inits*/
  double fx;
  double dfx[omp.D];
  double FI[omp.D*omp.D];
  double dFI[omp.D*omp.D*omp.D];
  double x[omp.D];
  for (i=0;i<omp.D;i++) x[i]=1;
  Posterior(x, &omp, &fx, dfx, FI, dFI );
  printf("x=ones(D); Posterior(x)=%f;\n",fx);
  PosteriorFI(x, &omp, FI);
  printf("[FI]\n");
  for (i=0;i<omp.D;i++){
    for (j=0;j<omp.D;j++) printf("%f\t",FI[i*omp.D+j]);
    printf("\n");
  }
  printf("[/FI]\n#\n");

  printf("[dFI]\n");
  for (k=0;k<omp.D;k++){
    for (j=0;j<omp.D;j++){
      for (i=0;i<omp.D;i++) printf("%f\t",dFI[k*omp.D*omp.D+j*omp.D+i]);
      printf("\n");
    }
    printf("\n");
  }
  printf("[/dFI]\n");

  /* print first sample, initial values in init_x */
  /* mcmc_print_sample(kernel, stdout); */

  size_t acc_c = 0;
  double acc_rate;
  size_t it;
  int acc;
	
  oFile=fopen(sample_file,"w");

  size_t BurnInSamples = cnf_options.sample_size/10;
  fprintf(stdout, "Burn-in phase.\n");
  /* Burn-in loop */
  for (it = 0; it < BurnInSamples; it++) {
    /* draw a sample using RMHMC */
    mcmc_sample(kernel, &acc);
    acc_c += acc;
    /* print sample */
    mcmc_print_sample(kernel, stdout);
		
    /* Addapt RMHMC parameters every 100 burn-in samples */
    if ( ((it + 1) % 100) == 0 ) {
      acc_rate = (double)acc_c / (double)100;
      fprintf(stdout, "Iteration: %li  Acceptance rate: %.2g.\n",it, acc_rate);
      /* print RMHM log and statistics */
      mcmc_print_stats(kernel, stdout);
      /* adapt RMHMC parameters */
      mcmc_adapt(kernel, acc_rate);
      acc_c = 0;
    }
  }
	
  fprintf(stdout, "\nBurn-in complete, sampling from the posterior.\n");
  /* Posterior loop */
  clock_t ct=clock();
  size_t Samples = cnf_options.sample_size;
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
      fprintf(stdout, "Iteration: %li  Acceptance rate: %.2g \n",it, acc_rate);
      mcmc_print_stats(kernel, stdout);
      acc_c = 0;
    }
  }
  ct=clock()-ct;
  printf("computation time spend sampling: %f s\n",((double) ct)/((double) CLOCKS_PER_SEC));

  fclose(oFile);
  /* clear memory */
  rmhmc_model_free(model);
  mcmc_free(kernel);
  ode_model_parameters_free(&omp);
  return 0 ;
}



/* Calculates the unormalised log-posterior fx, the gradient dfx, the Fisher information FI 
 * and the partial derivatives of the Fisher information matrix for ode models
 */
int Posterior(const double* x,  void* model_params, double* fx, double* dfx, double* FI, double* dFI ){
	
  ode_model_parameters* params =  (ode_model_parameters*) model_params;
  double prior_value=0;            // the prior distribution related sum of squares
  int i,j,F,T,D,U,status;
  gsl_matrix *Data,*sd;
  gsl_vector *fy;
  gsl_matrix *inv_cov;
  gsl_vector *prior_diff, *prior_mv;
  
  U=params->input_u->size2;
  D=params->D;
  prior_diff=params->prior_tmp_a;
  prior_mv=params->prior_tmp_b;

  fy=params->fy;
  Data=params->Data;
  sd=params->sdData;

  inv_cov=params->prior_inverse_cov;

  T=(params->t)->size;
  F=fy->size;
  // init FI and dFI
  fx[0]=0;

  gsl_sf_result exp_x;
  

  for (i=0;i<D;i++) {
    status=gsl_sf_exp_e(x[i],&exp_x);
    if (status==GSL_SUCCESS){
      if (exp_x.err/exp_x.val>0.1) printf("[exp] relative error larger than 10%%\n");
      gsl_vector_set(params->exp_x_u,i,exp_x.val);
    }
    else {
      printf("[gsl_exp returned %i] with x[%i]=%f\n",status,i,x[i]);
      exit(status);    
    }

  }
  if (FI!=NULL){
    for (i=0;i<D*D;i++) FI[i]=0;
  }
  if (dFI!=NULL){
    for (i=0;i<D*D*D;i++) dFI[i]=0;
  }
  //  printf("calculating Likelihood with:\n");
  //  gsl_printf("exp_x",params->exp_x_u,0);
  Likelihood(params, fx, dfx, FI, dFI);

  for (i=0;i<D;i++) gsl_vector_set(prior_diff,i,x[i]-gsl_vector_get(params->prior_mu,i));

  gsl_blas_dsymv(CblasUpper,-0.5,inv_cov,prior_diff,0.0,prior_mv);
  gsl_blas_ddot(prior_diff,prior_mv,&prior_value);
  //  printf("PosteriorFI: prior_value=%f\n",prior_value);
  fx[0]+=prior_value;

  // the prior related part of the gradient of fx:
  if (dfx!=NULL){
    for (i=0; i < D ; i++)
      dfx[i]+= gsl_vector_get(prior_mv,i);
  }
  if (FI!=NULL){
    for (i=0; i < D*D ; i++)
      FI[i] += params->prior_inverse_cov->data[i];
  }

  return 0;
}

/* Calculates the Fisher information matrix for a multivariate normal.
 */
int PosteriorFI(const double* x,  void* model_params, double* FI){
  double fx;
  Posterior(x,model_params,&fx,NULL,FI,NULL);
  return EXIT_SUCCESS;
}
