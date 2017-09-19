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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _GNU_SOURCE
#include <libgen.h>
#endif
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
#include "../mcmc/smmala_posterior.h"
#include "diagnosis_output.h"

// sampling actions 
#define SMPL_RESUME 1
#define SMPL_FRESH 0
#define SMPL_RESUME_TUNE 2

#define BUFSZ 1024
//#define BETA(rank,R) gsl_pow_4((double)(rank)/(double) ((R)-1))
//#define BETA(rank,R) (1.0/((double)(rank+1)))
#define BETA(rank,R) gsl_sf_exp(-0.25*((double) rank))

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
  int c,i=0;
  int warm_up=0; // sets the number of burn in points at command line
  char *cfilename=NULL;
  char lib_name[BUFSZ];
  ode_model_parameters omp;
  FILE *cnf;   // configuration file, with file name: cfilename
  char global_sample_filename_stem[BUFSZ]="sample.dat"; // filename basis
  char rank_sample_file[BUFSZ]; // filename for sample output
  //char *x_sample_file=NULL; // filename for sample output x(t,p)
  //char *y_sample_file=NULL; // filename for sample output y(t,p)
  FILE *oFile; // will be the file named «sample_file»
  FILE *rFile; // last sampled value will be written to this file
  char resume_filename[BUFSZ]="resume.double";
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
  int rank,R,DEST;
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
  
  char *dot;
  char *lib_base;
  lib_base=basename(lib_name);
  dot=strchr(lib_base,'.');
  dot[0]='\0';
  
  sprintf(resume_filename,"%s_resume_%04i.double",lib_base,rank);
  sprintf(rank_sample_file,"rank_%04i_of_%i_%s_%s",rank,R,lib_base,basename(cnf_options.output_file));
  cnf_options.output_file=rank_sample_file;
  
  printf("#\tlib_name: %s\nrank_sammple_file: %s\n#\tresume_filename: %s\n",lib_name,rank_sample_file,resume_filename);
  
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
  //  if(rank==0)  gsl_printf("default parameters",&(p_view.vector),GSL_IS_DOUBLE | GSL_IS_VECTOR);
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
  
  /* if (rank==0){ */
  /*   gsl_printf("data",omp.Data,GSL_IS_DOUBLE | GSL_IS_MATRIX); */
  /*   gsl_printf("sd data",omp.sdData,GSL_IS_DOUBLE | GSL_IS_MATRIX); */
  /*   //gsl_vector_fprintf(stdout,omp.E[0]->input_u,"%g"); fflush(stdout); */
  /*   //printf("Experiment %i\n",c); */
  /*   for (c=0;c<omp.size->C;c++) gsl_printf("u",omp.E[c]->input_u,GSL_IS_DOUBLE | GSL_IS_VECTOR); */
  /*   for (c=0;c<omp.size->C;c++) gsl_printf("t",omp.E[c]->t,GSL_IS_DOUBLE | GSL_IS_VECTOR); */
  /* } */
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
  smmala_model* model = smmala_model_alloc(LogPosterior, NULL, &omp);
  
  /* initial parameter values */
  D=omp.size->D;
  double init_x[D];
  /* allocate a new RMHMC MCMC kernel */
  printf("# allocating memory for a new SMMALA MCMC kernel.\n");
  double beta=BETA(rank,R);
  mcmc_kernel* kernel = smmala_kernel_alloc(beta,D,
					    cnf_options.initial_stepsize,
					    model,
					    seed,
					    cnf_options.target_acceptance);
  /* initialise MCMC */
  // if (rank==0) gsl_printf("prior mean",omp.prior_mu,GSL_IS_DOUBLE | GSL_IS_VECTOR);  
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
  /* if (rank==0){ */
  /*   gsl_printf("Normalised Data",omp.Data,GSL_IS_DOUBLE | GSL_IS_MATRIX); */
  /*   gsl_printf("standard deviation",omp.sdData,GSL_IS_DOUBLE | GSL_IS_MATRIX); //exit(0); */
  /*   fflush(stdout); */
  /* } */
  printf("# test evaluation of Posterior function done:\n");
  printf("# θ=θ₀; Posterior(θ|D)=%+g;\n# where θ₀:",kernel->fx[0]); 
  for (i=0;i<D; i++) printf(" %+g ",kernel->x[i]); printf("\n");

  void *buffer=(void *) smmala_comm_buffer_alloc(D);

  /* print first sample, initial values in init_x */
  mcmc_print_sample(kernel, stdout);
  ode_solver_print_stats(solver, stdout);
  fflush(stdout);
  fflush(stderr);
  
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
  printf("# Performing Burn-In with step-size (%g) tuning: %lu iterations\n",cnf_options.initial_stepsize,BurnInSamples);
  int master=0;
  int swaps=0;
  /* Burn In Loop and Tuning*/
  for (it = 0; it < BurnInSamples; it++) {
    mcmc_sample(kernel, &acc);
    acc_c += acc;
    if ((it+1)%3==0){
      master=(it%2==rank%2); // if iterator is even, then even procs are master
      if (master){
	DEST=(rank+1)%R; // this process is the master process for swap decisions
      } else {
	DEST=(R+rank-1)%R; // this process has to accept swap decisions from DEST
      }
      mcmc_exchange_information(kernel,DEST,buffer);
      swaps+=mcmc_swap_chains(kernel,master,rank,DEST,buffer);
    }
    //mcmc_print_sample(kernel, stdout);
    if ( ((it + 1) % 100) == 0 ) {
      acc_rate = (double) acc_c / 100.0;
      fprintf(stdout, "# [rank %i/%i; β=%5f] (it %li) acc. rate: %3.2g; %2i %% swaps\t",rank,R,beta,it,acc_rate,swaps);
      mcmc_print_stats(kernel, stdout);
      mcmc_adapt(kernel, acc_rate);
      acc_c = 0;
      swaps=0;
    }
  }
  
  fprintf(stdout, "\n# Burn-in complete, sampling from the posterior.\n");
  
  
  /* full Posterior loop */
  size_t Samples = cnf_options.sample_size;
  clock_t ct=clock();
  
  for (it = 0; it < Samples; it++) {
    /* draw a sample using RMHMC */
    mcmc_sample(kernel, &acc);
    acc_c += acc;
    master=(it%2==rank%2);
    if (master){
      DEST=(rank+1)%R; // this process is the master process for swap decisions
    } else {
      DEST=(R+rank-1)%R; // this process has to accept swap decisions from DEST
    }
    //their_beta=BETA(DEST,R); // the orther proc's 
    mcmc_exchange_information(kernel,DEST,buffer);
    mcmc_swap_chains(kernel,master,rank,DEST,buffer);
		
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



