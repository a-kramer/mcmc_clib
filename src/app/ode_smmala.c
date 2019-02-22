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
#include <assert.h>
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
#include "flatten.h"
#include "read_data_hdf5.h"
#include "../mcmc/smmala.h"
#include "../ode/ode_model.h"
#include "../mcmc/smmala_posterior.h"
#include "diagnosis_output.h"
#include "hdf5.h"
#include "hdf5_hl.h"

#include "../mcmc/model_parameters_smmala.h"
// define target block types
#define INTEGER_BLOCK 1
#define DOUBLE_BLOCK 2

// data field ids, should be consecutive atm., because they are sometimes looped over.
#define i_time 0
#define i_reference_input 1
#define i_reference_data 2
#define i_sd_reference_data 3
#define i_input 4
#define i_data 5
#define i_sd_data 6
#define i_prior_mu 7
#define i_prior_icov 8
#define i_initial_conditions 9
#define i_ref_initial_conditions 10
#define i_norm_f 11
#define i_norm_t 12
#define NumberOfFields 13

typedef struct {
  char *library_file;
  char *output_file;
  double target_acceptance;
  double initial_stepsize;
  long sample_size;
  double t0;
} main_options;  // these are user supplied options to the program


#define CHUNK 100
// sampling actions 
#define SMPL_RESUME 1
#define SMPL_FRESH 0
#define SMPL_RESUME_TUNE 2

#define BUFSZ 2048
//#define BETA(rank,R) gsl_pow_4((double)(rank)/(double) ((R)-1))
//#define BETA(rank,R) (1.0/((double)(rank+1)))
#define BETA(rank,R) gsl_sf_exp(-gamma*((double) rank))

/* Auxiliary structure with working storage and aditional parameters for
 * a multivariate normal model with known covariance matrix and zero mean.
 *
typedef struct {
	int D
	double* Variance;
	double* Precision;
	double* tmpVec;
	char init;
} mvNormParams;
 */

int print_help(){
  printf("Usage:\n");
  printf("-a $ACCEPTANCE_RATE\n");
  printf("\t\t\tTarget acceptance value (all markov chains will be tuned for this acceptance).\n\n");
  printf("-d, --hdf5 ./data.h5\n");
  printf("\t\t\tdata.h5 contains the data points and the conditions of measurement in hdf5 format. This and .cfg files are mutually exclusive. A suitable h5 file is produced by the hdf5_import program.\n\n");
  printf("-g $G\n");
  printf("\t\t\tThis will define how the inverse MCMC temperatures β are chosen: β = exp(-G*MPI_Rank).\n\n");
  printf("-i $STEP_SIZE\n");
  printf("\t\t\tThe initial step size of each markov chain, this will usually be tuned to get the desired acceptance rate $A (-a $A).\n\n");
  printf("-l ./ode_model.so\n");
  printf("\t\t\tode_model.so is a shared library containing the CVODE functions of the model.\n\n");
  printf("-o ./output_file.h5\n");
  printf("\t\t\tFilename for hdf5 output. This file will contain the log-parameter sample and log-posterior values. The samples will have attributes that reflect the markov chain setup.\n\n");
  printf("-p, --prior-start\n");
  printf("\t\t\tStart the markov chain at the center of the prior.\n");
  printf("-r, --resume\n");
  printf("\t\t\tResume from last sampled MCMC point. Only the last MCMC position is read from the file named «resume.double». Everything else about the problem can be changed.\n");
  printf("-s $N\n");
  printf("\t\t\t$N sample size. default N=10.\n\n");
  printf("-t,--init-at-t $T_INITIAL\n");
  printf("\t\t\tSpecifies the initial time «t0» of the model integration [initial value problem for the ordinary differential equation in x; x(t0)=x0]\n\n");
  //printf("-b\n");
  //printf("\t\t\tOutput mode: binary.\n\n");
  printf("--seed $SEED\n");
  printf("\t\tSet the gsl pseudo random number generator seed to $SEED. (perhaps --seed $RANDOM)\n\n");
  MPI_Finalize();
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
  char global_sample_filename_stem[BUFSZ]="Sample.h5"; // filename basis
  char rank_sample_file[BUFSZ]; // filename for sample output
  //char *x_sample_file=NULL; // filename for sample output x(t,p)
  //char *y_sample_file=NULL; // filename for sample output y(t,p)
  FILE *oFile; // will be the file named «sample_file»
  FILE *rFile; // last sampled value will be written to this file
  char resume_filename[BUFSZ]="resume.double";
  //int output_is_binary=0;
  double seed = 1;
  double gamma=0.25;
  double t0=NAN;
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
  herr_t status;
  MPI_Init(&argc,&argv);
  int rank,R,DEST;
  MPI_Comm_size(MPI_COMM_WORLD,&R);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  char *h5file=NULL;
  for (i=0;i<argc;i++){
    if (strcmp(argv[i],"-p")==0 || strcmp(argv[i],"--prior-start")==0) {
      start_from_prior=1;
    } else if (strcmp(argv[i],"-d")==0 || strcmp(argv[i],"--hdf5")==0) {
      h5file=argv[i+1];
    } else if (strcmp(argv[i],"-t")==0 || strcmp(argv[i],"--init-at-t")==0) {
      t0=strtod(argv[i+1],NULL);
      printf("[main] t0=%f",t0);
    } else if (strcmp(argv[i],"-w")==0 || strcmp(argv[i],"--warm-up")==0) warm_up=strtol(argv[i+1],NULL,10);
    else if (strcmp(argv[i],"--resume")==0 || strcmp(argv[i],"-r")==0) sampling_action=SMPL_RESUME;
    else if (strcmp(argv[i],"--sens-approx")==0) sensitivity_approximation=1;
    else if (strcmp(argv[i],"-l")==0) strcpy(cnf_options.library_file,argv[i+1]);
    //    else if (strcmp(argv[i],"-n")==0) Tuning=0;
    else if (strcmp(argv[i],"-s")==0) cnf_options.sample_size=strtol(argv[i+1],NULL,0);
    else if (strcmp(argv[i],"-o")==0) strncpy(cnf_options.output_file,argv[i+1],BUFSZ);
    //    else if (strcmp(argv[i],"-b")==0) output_is_binary=1;
    else if (strcmp(argv[i],"-a")==0) cnf_options.target_acceptance=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-i")==0) cnf_options.initial_stepsize=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-g")==0) gamma=strtod(argv[i+1],NULL);
    
    else if (strcmp(argv[i],"--seed")==0) seed=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--help")==0) {print_help(); MPI_Abort(MPI_COMM_WORLD,0);}
    
    //else printf("unknown option {%s}.\n",argv[i]);
  }
  
  seed=seed*137+13*rank;
  
  /* load model */
  ode_model *odeModel = ode_model_loadFromFile(lib_name);  /* alloc */
  if (odeModel == NULL) {
    fprintf(stderr, "# [main] Library %s could not be loaded.\n",lib_name);
    exit(1);
  } else printf( "# [main] Library %s loaded.\n",lib_name);
  
  char *dot;
  char *lib_base;
  lib_base=basename(lib_name);
  dot=strchr(lib_base,'.');
  dot[0]='\0';
  sprintf(resume_filename,"%s_resume_%02i.double",lib_base,rank);
  sprintf(rank_sample_file,"mcmc_rank_%02i_of_%i_%s_%s",rank,R,lib_base,basename(cnf_options.output_file));
  cnf_options.output_file=rank_sample_file;
  
  //printf("#\tlib_name: %s\nrank_sammple_file: %s\n#\tresume_filename: %s\n",lib_name,rank_sample_file,resume_filename);
  
  ode_solver* solver = ode_solver_alloc(odeModel); /* alloc */
  if (solver == NULL) {
    fprintf(stderr, "# [main] Solver %s could not be created.\n",lib_base);
    ode_model_free(odeModel);
    exit(1);
  } else printf("# [main] Solver for %s created.\n",lib_base);
  
  if (sensitivity_approximation){
    ode_solver_disable_sens(solver);
    odeModel->vf_sens=NULL; // make sens function unavailable; ode_model_has_sens(model) will return «False»;
  }
  /* init solver */
  realtype solver_param[3] = {ODE_SOLVER_ABS_ERR, ODE_SOLVER_REL_ERR, 0};

  char **x_name, **p_name, **f_name;
  x_name=ode_model_get_var_names(odeModel);
  p_name=ode_model_get_param_names(odeModel);
  f_name=ode_model_get_func_names(odeModel);
  
  /* define local variables for parameters and inital conditions */
  int N = ode_model_getN(odeModel);
  int P = ode_model_getP(odeModel);
  int F = ode_model_getF(odeModel);
  omp.size=(problem_size*) malloc(sizeof(problem_size));
  omp.size->N=N;
  omp.size->P=P;
  omp.size->F=F;
  omp.t0=t0;
  double y[N];
  gsl_vector_view y_view=gsl_vector_view_array(y,N);
  
  double p[P];
  gsl_vector_view p_view=gsl_vector_view_array(p,P);
  ode_model_get_default_params(odeModel, p, P);
  if(rank==0)  gsl_printf("default parameters",&(p_view.vector),GSL_IS_DOUBLE | GSL_IS_VECTOR);
  //  double t_old = solver_param[2];
  
  omp.solver=solver;

  if (h5file!=NULL){
    printf("# [main] reading hdf5 file, loading data..."); fflush(stdout);
    read_data(h5file,&omp);
    printf("done.\n"); fflush(stdout);
  } else {
    printf("# [main] no data provided (-c or -d option), exiting.\n");
    MPI_Abort(MPI_COMM_WORLD,-1);
  }
  
  cnf_options.initial_stepsize=fabs(cnf_options.initial_stepsize);
  cnf_options.target_acceptance=fabs(cnf_options.target_acceptance);
  cnf_options.sample_size=fabs(cnf_options.sample_size);
  ode_model_get_initial_conditions(odeModel, y, N);
  /*
   * if (omp.normalisation_type==DATA_IS_ABSOLUTE){
   *  ode_model_get_initial_conditions(odeModel, y, N);  
   * } else {
   *  gsl_vector_set_all(&(y_view.vector),1.0);
   * }
   */
  int C=omp.size->C;
  
  /* if (rank==0){ */
  /*   for (c=0;c<C;c++){ */
  /*     printf("[main] Experiment %i:\n",c); */
  /*     gsl_printf("data",omp.E[c]->data_block,GSL_IS_DOUBLE | GSL_IS_MATRIX); */
  /*     gsl_printf("standard deviation",omp.E[c]->sd_data_block,GSL_IS_DOUBLE | GSL_IS_MATRIX); */
  /*     gsl_printf("u",omp.E[c]->input_u,GSL_IS_DOUBLE | GSL_IS_VECTOR); */
  /*     gsl_printf("t",omp.E[c]->t,GSL_IS_DOUBLE | GSL_IS_VECTOR); */
  /*   } */
  /* } */
  // unspecified initial conditions
  
  int NNE=0; // number of normalising experiments
  int LE=0;  
  for (i=0;i<C;i++){
    if (omp.E[i]->init_y==NULL){
      omp.E[i]->init_y=gsl_vector_alloc(N);
      gsl_vector_memcpy(omp.E[i]->init_y,&(y_view.vector));
    }
    if (omp.E[i]->lflag==0) NNE++;
  }
  
  LE=C-NNE;
  if (rank==0){
    printf("# [main] There are %i experiments",C);
    if (NNE>0){
      printf(", %i of which ",NNE);
      if (NNE==1) printf("is");
      else printf("are");
      printf(" used only for the normalisation of the %i experiments that explicitly contribute to the LogLikelihood(NormalisedData[1:%i]|θ).\n",LE,LE); 	
    }  else printf(".\n");
  }    
  printf("# [main] init ivp: t0=%g\n",omp.t0); fflush(stdout);
  ode_solver_init(solver, omp.t0, y, N, p, P);
  printf("# [main] solver initialised.\n"); fflush(stdout);
  ode_solver_setErrTol(solver, solver_param[1], &solver_param[0], 1);
  if (ode_model_has_sens(odeModel)) {
    ode_solver_init_sens(solver, omp.E[0]->yS0->data, P, N);
    printf("# [main] sensitivity analysis initiated.\n");
  }
  
  //printf("# [main] allocating memory for the SMMALA model.\n");
  fflush(stdout);
  smmala_model* model = smmala_model_alloc(LogPosterior, NULL, &omp);
  if (!model){
    fprintf(stderr,"smmala_model could not be allocated.");
    MPI_Abort(MPI_COMM_WORLD,-1);
  }
  /* initial parameter values */
  D=omp.size->D;
  double init_x[D];
  /* allocate a new RMHMC MCMC kernel */
  //printf("# [main] allocating memory for a new SMMALA MCMC kernel.\n");
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
      fprintf(stderr,"[rank %i] Could not open resume file. Starting from: ",rank);
      for (i=0;i<D;i++) init_x[i]=gsl_vector_get(omp.prior->mu,i);
    } else {
      resume_count=fread(init_x, sizeof(double), D, rFile);
      fclose(rFile);
      if (resume_count!=D) {
	fprintf(stderr,"[rank %i] Reading from resume file returned a wrong number of values.",rank);
	exit(-1);
      }
    }
  } else if (start_from_prior==1){     
    if (rank==0) printf("# [main] setting initial mcmc vector to prior mean.\n");
    for (i=0;i<D;i++) init_x[i]=gsl_vector_get(omp.prior->mu,i);
  } else {
    if (rank==0) printf("# [main] setting mcmc initial value to log(default parameters)\n");
    for (i=0;i<D;i++) init_x[i]=gsl_sf_log(p[i]);
  }
  //display_prior_information(omp.prior);
  if (rank==0){
    printf("# [main] initializing MCMC.\n");
    printf("# [main] init_x:");
    for (i=0;i<D;i++) printf(" %g ",init_x[i]);
    printf("\n");
  }
  mcmc_init(kernel, init_x);
  printf("# [main] rank %i init complete .\n",rank);
   
  //gsl_printf("mu",omp.prior->mu,GSL_IS_DOUBLE|GSL_IS_VECTOR);
  if (rank==0){
    printf("# [main] test evaluation of Posterior function done:\n");
    printf("# \tθ=θ₀; LogPosterior(θ|D)=%+g;\n# where θ₀:",kernel->fx[0]); 
    for (i=0;i<D;i++) printf(" %+g ",kernel->x[i]); printf("\n");
    printf("# [main] LogLikelihood(D|θ):");
    printf("%+g\tLogPrior(θ)=%+g.\n",kernel->fx[1],kernel->fx[2]);    
  }
  ode_solver_print_stats(solver, stdout);
  fflush(stdout);
  fflush(stderr);
  MPI_Abort(MPI_COMM_WORLD,0);
  
  void *buffer=(void *) smmala_comm_buffer_alloc(D);
  size_t acc_c = 0;
  double acc_rate;
  size_t it;
  int acc;
  size_t Samples = cnf_options.sample_size;  
  //oFile=fopen(cnf_options.output_file,"w");
  hsize_t size[2];
  hsize_t chunk_size[2];
  hid_t file_id = H5Fcreate(cnf_options.output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  assert(file_id>0);

  char *x_names=flatten(x_name, (size_t) N, "; ");
  char *p_names=flatten(p_name, (size_t) P, "; ");
  char *f_names=flatten(f_name, (size_t) F, "; ");

  herr_t NameWriteError=0;
  NameWriteError&=H5LTmake_dataset_string(file_id,"StateVariableNames",x_names);
  NameWriteError&=H5LTmake_dataset_string(file_id,"ParameterNames",p_names);
  NameWriteError&=H5LTmake_dataset_string(file_id,"OutputFunctionNames",f_names);
  if (NameWriteError==0) printf("# [main] all names have been written to file «%s».\n",cnf_options.output_file);
  free(x_names);
  free(p_names);
  free(f_names);
  
  hid_t para_property_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t post_property_id = H5Pcreate(H5P_DATASET_CREATE);
  
  /* here we set a «chunk size», which will coincide with the hyperslabs we select to write output*/
  // parameter sample chunk size:
  chunk_size[0]=CHUNK;
  chunk_size[1]=D;
  H5Pset_chunk(para_property_id, 2, chunk_size);
  hid_t para_chunk_id=H5Screate_simple(2, chunk_size, NULL); // a dataspace to write chunks/hyperslabs
  
  // posterior probability distribution chunk size:
  chunk_size[0]=CHUNK;
  chunk_size[1]=1;
  H5Pset_chunk(post_property_id, 2, chunk_size);
  hid_t post_chunk_id=H5Screate_simple(2, chunk_size, NULL);
  // hyperslab selection
  hsize_t offset[2]={0,0}, stride[2]={1,1}, count[2]={1,1}, block[2];
  block[0]=CHUNK;
  block[1]=D;

  // hdf5 file setup
  size[0]=Samples;
  size[1]=D;
  hid_t para_dataspace_id=H5Screate_simple(2, size, NULL);
  size[0]=Samples;
  size[1]=1;
  hid_t post_dataspace_id=H5Screate_simple(2, size, NULL);
    assert(post_dataspace_id>0 && para_dataspace_id>0);
  
  hid_t parameter_set_id = H5Dcreate2(file_id, "LogParameters", H5T_NATIVE_DOUBLE, para_dataspace_id, H5P_DEFAULT, para_property_id, H5P_DEFAULT);
  hid_t posterior_set_id = H5Dcreate2(file_id, "LogPosterior", H5T_NATIVE_DOUBLE, post_dataspace_id, H5P_DEFAULT, post_property_id, H5P_DEFAULT);
  assert(parameter_set_id>0 && posterior_set_id>0);

  // Burn In Loop set up
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
      if (R>1){
	mcmc_exchange_information(kernel,DEST,buffer);
	swaps+=mcmc_swap_chains(kernel,master,rank,DEST,buffer);
      }
    }
    //mcmc_print_sample(kernel, stdout);
    if ( ((it + 1) % CHUNK) == 0 ) {
      acc_rate = ((double) acc_c) / ((double) CHUNK);
      fprintf(stdout, "# [rank %i/%i; β=%5f] (it %4li) acc. rate: %3.2f; %3i %% swap success\t",rank,R,beta,it,acc_rate,swaps);
      mcmc_print_stats(kernel, stdout);
      mcmc_adapt(kernel, acc_rate);
      acc_c = 0;
      swaps=0;
    }
  }  
  fprintf(stdout, "\n# Burn-in complete, sampling from the posterior.\n");
  /* full Posterior loop */
  clock_t ct=clock();
  gsl_matrix *log_para_chunk;
  gsl_vector *log_post_chunk;
  log_para_chunk=gsl_matrix_alloc(CHUNK,D);
  log_post_chunk=gsl_vector_alloc(CHUNK);
  gsl_vector_view current;
  gsl_vector_view x_state;
  swaps=0;
  acc_c = 0;
  for (it = 0; it < Samples; it++) {
    mcmc_sample(kernel, &acc);
    acc_c += acc;
    master=(it%2==rank%2);
    if (master){
      DEST=(rank+1)%R; // this process is the master process for swap decisions
    } else {
      DEST=(R+rank-1)%R; // this process has to accept swap decisions from DEST
    }
    //their_beta=BETA(DEST,R); // the other proc's
    if (R>1){
      mcmc_exchange_information(kernel,DEST,buffer);
      swaps+=mcmc_swap_chains(kernel,master,rank,DEST,buffer);
    }
    /* save sampled point for writing */
    current=gsl_matrix_row(log_para_chunk,it%CHUNK);
    x_state=gsl_vector_view_array(kernel->x,D);
    gsl_vector_memcpy(&(current.vector),&(x_state.vector));
    gsl_vector_set(log_post_chunk,it%CHUNK,kernel->fx[0]);
    /* print sample log and statistics every 100 samples */
    if ( ((it + 1) % CHUNK) == 0 ) {
      acc_rate = ((double) acc_c) / ((double) CHUNK);
      fprintf(stdout, "# [rank %i/%i; β=%5f; %3li%% done] (it %5li) acc. rate: %3.2f; %3i %% swap success\t",rank,R,beta,it/Samples,it,acc_rate,swaps);
      mcmc_print_stats(kernel, stdout);
      acc_c = 0;

      // print chunk to hdf5 file
      block[1]=D;
      status = H5Sselect_hyperslab(para_dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
      H5Dwrite(parameter_set_id, H5T_NATIVE_DOUBLE, para_chunk_id, para_dataspace_id, H5P_DEFAULT, log_para_chunk->data);
      block[1]=1;
      status = H5Sselect_hyperslab(post_dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
      H5Dwrite(posterior_set_id, H5T_NATIVE_DOUBLE, post_chunk_id, post_dataspace_id, H5P_DEFAULT, log_post_chunk->data);
      offset[0]=it+1;
      swaps=0;
    }
  }
  // write remaining data to the output hdf5 file
  int Rest=Samples % CHUNK;
  printf("[main] last iteration done %i points remain to write.\n",Rest);
  if (Rest > 0){
    chunk_size[0]=Rest;
    chunk_size[1]=D;
    para_chunk_id=H5Screate_simple(2, chunk_size, NULL);    
    
    chunk_size[0]=Rest;
    chunk_size[1]=1;
    post_chunk_id=H5Screate_simple(2, chunk_size, NULL);
    
    printf("[main] writing the remaining %i sampled parametrisations to file.\n",Rest);
    block[0]=Rest;
    block[1]=D;
    printf("[main] offset: %lli×%lli; block: %lli×%lli; stride: %lli×%lli; count: %lli×%lli.\n",offset[0],offset[1],block[0],block[1],stride[0],stride[1],count[0],count[1]);
    status = H5Sselect_hyperslab(para_dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
    H5Dwrite(parameter_set_id, H5T_NATIVE_DOUBLE, para_chunk_id, para_dataspace_id, H5P_DEFAULT, log_para_chunk->data);
    block[1]=1;
    printf("[main] writing their %i log-posterior values to file.\n",Rest);
    printf("[main] offset: %lli×%lli; block: %lli×%lli; stride: %lli×%lli; count: %lli×%lli.\n",offset[0],offset[1],block[0],block[1],stride[0],stride[1],count[0],count[1]);
    status &= H5Sselect_hyperslab(post_dataspace_id, H5S_SELECT_SET, offset, stride, count, block);
    H5Dwrite(posterior_set_id, H5T_NATIVE_DOUBLE, post_chunk_id, post_dataspace_id, H5P_DEFAULT, log_post_chunk->data);
    assert(status>=0);
  }

  // annotate written sample with all necessary information
  printf("[main] writing some annotation about the sampled points as hdf5 attributes.\n");
  status&=H5LTset_attribute_double(file_id, "LogParameters", "seed", &seed, 1);
  status&=H5LTset_attribute_int(file_id, "LogParameters", "MPI_RANK", &rank, 1);
  status&=H5LTset_attribute_ulong(file_id, "LogParameters", "SampleSize", &Samples, 1);
  status&=H5LTset_attribute_ulong(file_id, "LogParameters", "BurnIn", &BurnInSamples, 1);
  status&=H5LTset_attribute_double(file_id, "LogParameters", "Temperature", &beta, 1);
  
  status&=H5LTset_attribute_string(file_id, "LogParameters", "ModelLibrary", lib_base);
  status&=H5LTset_attribute_string(file_id, "LogParameters", "DataFrom", cfilename);
  ct=clock()-ct;
  printf("# computation time spend sampling: %f s\n",((double) ct)/((double) CLOCKS_PER_SEC));
  
  rFile=fopen(resume_filename,"w");
  mcmc_write_sample(kernel, rFile);
  fclose(rFile);

  H5Dclose(posterior_set_id);
  H5Dclose(parameter_set_id);
  H5Sclose(para_dataspace_id);
  H5Sclose(post_dataspace_id);
  H5Pclose(para_property_id);
  H5Pclose(post_property_id);
  H5Fclose(file_id);

  /* clear memory */
  smmala_model_free(model);
  mcmc_free(kernel);
  ode_model_parameters_free(&omp);
  MPI_Finalize();
  return EXIT_SUCCESS;
}



