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

#define yes 1
#define no 0
#define CHUNK 100
// sampling actions 
#define SMPL_RESUME 1
#define SMPL_FRESH 0
#define SMPL_RESUME_TUNE 2

#define BUFSZ 2048
//#define BETA(rank,R) gsl_pow_4((double)(R-rank)/(double) (R))
//#define BETA(rank,R) (1.0/((double)(rank+1)))
//#define BETA(rank,R) gsl_sf_exp(-gamma*((double) rank))



typedef struct {
  char *library_file;
  char *output_file;
  char *resume_file;
  double target_acceptance;
  double initial_stepsize;
  double initial_stepsize_rank_factor; // step_size = pow(rank_factor,rank) * initial_step_size;
  long sample_size;
  double abs_tol;
  double rel_tol;
  double t0;
} main_options;  // these are user supplied options to the program

typedef struct {  
  hid_t file_id;
  hid_t para_property_id;
  hid_t post_property_id;
  hid_t para_chunk_id;
  hid_t post_chunk_id;
  hid_t para_dataspace_id;
  hid_t post_dataspace_id;
  hid_t parameter_set_id;
  hid_t posterior_set_id;
  hsize_t *size;
  hsize_t *chunk_size;
  hsize_t *offset;
  hsize_t *stride;
  hsize_t *count;
  hsize_t *block;
} hdf5block_t;

int h5block_close(hdf5block_t *h5block){
  H5Dclose(h5block->posterior_set_id);
  H5Dclose(h5block->parameter_set_id);
  H5Sclose(h5block->para_dataspace_id);
  H5Sclose(h5block->post_dataspace_id);
  H5Pclose(h5block->para_property_id);
  H5Pclose(h5block->post_property_id);
  H5Fclose(h5block->file_id);
  free(h5block->size);
  free(h5block->chunk_size);
  free(h5block->offset);
  free(h5block->stride);
  free(h5block->count);
  free(h5block->block);
  free(h5block);
  return EXIT_SUCCESS;
}

hdf5block_t* h5block_init(char *output_file, ode_model_parameters *omp, size_t Samples, const char **x_name, const char **p_name, const char **f_name){
  hsize_t *size=malloc(sizeof(hsize_t)*2);
  hsize_t *chunk_size=malloc(sizeof(hsize_t)*2);
  hid_t file_id = H5Fcreate(output_file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  assert(file_id>0);
  int D=get_number_of_MCMC_variables(omp);
  int N=get_number_of_state_variables(omp);
  int P=get_number_of_model_parameters(omp);
  int F=get_number_of_model_outputs(omp);
  
  char *x_names=flatten(x_name, (size_t) N, "; ");
  char *p_names=flatten(p_name, (size_t) P, "; ");
  char *f_names=flatten(f_name, (size_t) F, "; ");

  herr_t NameWriteError=0;
  NameWriteError&=H5LTmake_dataset_string(file_id,"StateVariableNames",x_names);
  NameWriteError&=H5LTmake_dataset_string(file_id,"ParameterNames",p_names);
  NameWriteError&=H5LTmake_dataset_string(file_id,"OutputFunctionNames",f_names);
  if (NameWriteError){
    fprintf(stderr,"[h5block_init] writing (x,p,f)-names into hdf5 file failed.");
  }/* else {
    printf("# [main] all names have been written to file «%s».\n",output_file);
    }*/
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
  hsize_t *offset, *stride, *count, *block;
  offset=malloc(sizeof(hsize_t)*2);
  stride=malloc(sizeof(hsize_t)*2);
  count=malloc(sizeof(hsize_t)*2);
  block=malloc(sizeof(hsize_t)*2);  
  // hsize_t offset[2]={0,0}, stride[2]={1,1}, count[2]={1,1}, block[2];
  int i;
  for (i=0;i<2;i++) offset[i]=0;
  for (i=0;i<2;i++) stride[i]=1;
  for (i=0;i<2;i++)  count[i]=1;  
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
  //copy ids to struct
  hdf5block_t *h5block;
  h5block = malloc(sizeof(hdf5block_t));
  h5block->file_id = file_id;
  h5block->para_property_id = para_property_id;
  h5block->post_property_id = post_property_id;
  h5block->para_chunk_id = para_chunk_id;
  h5block->post_chunk_id = post_chunk_id;
  h5block->para_dataspace_id = para_dataspace_id;
  h5block->post_dataspace_id = post_dataspace_id;
  h5block->parameter_set_id = parameter_set_id;
  h5block->posterior_set_id = posterior_set_id;
  h5block->size = size;
  h5block->chunk_size = chunk_size;
  h5block->offset = offset;
  h5block->stride = stride;
  h5block->count = count;
  h5block->block = block;
  return h5block;
}

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

double assign_beta(int rank, int R, int gamma){
  double x=(double)(R-rank)/(double) R;
  double b=-1;
  assert(gamma>=1);
  b=gsl_pow_int(x, gamma);
  assert(b>=0 && b<=1.0);
  return b;
}

void print_help(){
  printf("Usage ($SOMETHING are values you choose, written as bash variables):\n");
  printf("-a $ACCEPTANCE_RATE\n");
  printf("\t\t\tTarget acceptance value (all markov chains will be tuned for this acceptance).\n\n");
  printf("-d, --hdf5 ./data.h5\n");
  printf("\t\t\tdata.h5 is a file that contains the data points and the conditions of measurement in hdf5 format. A suitable h5 file is produced by the hdf5_import program bundled with ode_smmala.\n\n");
  printf("-g $G\n");
  printf("\t\t\tThis will define how the inverse MCMC temperatures β are chosen: β = (1-rank/R)^G, where R is MPI_Comm_Size.\n\n");
  printf("-i $STEP_SIZE\n");
  printf("\t\t\tThe initial step size of each markov chain, this will usually be tuned to get the desired acceptance rate $A (-a $A).\n\n");
  printf("-l ./ode_model.so\n");
  printf("\t\t\tode_model.so is a shared library containing the CVODE functions of the model.\n\n");
  printf("-o ./output_file.h5\n");
  printf("\t\t\tFilename for hdf5 output. This file will contain the log-parameter sample and log-posterior values. The samples will have attributes that reflect the markov chain setup.\n\n");
  printf("-p, --prior-start\n");
  printf("\t\t\tStart the markov chain at the center of the prior. Otherwise it will be started from the DefaultParameters in the vfgen file.\n");
  printf("-r, --resume\n");
  printf("\t\t\tResume from last sampled MCMC point. Only the last MCMC position is read from the resume file. Everything else about the problem can be changed.\n");
  printf("-s $N\n");
  printf("\t\t\t$N sample size. default N=10.\n\n");
  printf("-t,--init-at-t $T_INITIAL\n");
  printf("\t\t\tSpecifies the initial time «t0» of the model integration [initial value problem for the ordinary differential equation in x; x(t0)=x0]\n\n");
  printf("--seed $SEED\n");
  printf("\t\tSet the gsl pseudo random number generator seed to $SEED. (perhaps --seed $RANDOM)\n\n");
  exit(EXIT_SUCCESS);
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
  for (j=0;j<width-2;j++) printf("-");
  printf("|\n");
  printf("%+4.4g",min);
  for (j=0;j<width-8;j++) printf(" ");
  printf("%+4.4g\n",max);  
}

void display_chunk_properties(hdf5block_t *h5block){
  hsize_t *offset=h5block->offset;
  hsize_t *block=h5block->block;
  hsize_t *count=h5block->count;
  hsize_t *stride=h5block->stride;
  printf("#[main]");
  printf("# offset: %lli×%lli;\n",offset[0],offset[1]);
  printf("#  block: %lli×%lli;\n",block[0],block[1]);
  printf("#  stride: %lli×%lli;\n",stride[0],stride[1]);
  printf("#  count: %lli×%lli.\n",count[0],count[1]);
}

herr_t h5write_current_chunk(hdf5block_t *h5block, gsl_matrix *log_para_chunk, gsl_vector *log_post_chunk){
  herr_t status;
  assert(log_para_chunk);
  assert(log_post_chunk);
  int D=log_para_chunk->size2;
  
  h5block->block[0]=CHUNK;
  h5block->block[1]=D;
  status = H5Sselect_hyperslab(h5block->para_dataspace_id, H5S_SELECT_SET, h5block->offset, h5block->stride, h5block->count, h5block->block);
  H5Dwrite(h5block->parameter_set_id, H5T_NATIVE_DOUBLE, h5block->para_chunk_id, h5block->para_dataspace_id, H5P_DEFAULT, log_para_chunk->data);

  h5block->block[1]=1;
  status = H5Sselect_hyperslab(h5block->post_dataspace_id, H5S_SELECT_SET, h5block->offset, h5block->stride, h5block->count, h5block->block);
  H5Dwrite(h5block->posterior_set_id, H5T_NATIVE_DOUBLE, h5block->post_chunk_id, h5block->post_dataspace_id, H5P_DEFAULT, log_post_chunk->data);
  return status;
}

int burn_in_foreach(int rank, int R, size_t BurnInSampleSize, ode_model_parameters *omp, mcmc_kernel *kernel, void *buffer){
  int master=0;
  int swaps=0;
  int acc=0, acc_c=0;
  double acc_rate=0.0;
  size_t it;
  double beta=mcmc_get_beta(kernel);
  int DEST;

  /* Burn In Loop and Tuning*/
  for (it = 0; it < BurnInSampleSize; it++) {
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
      fprintf(stdout, "# [rank % 2i/% 2i; β=%5f] (it %4li)\tacc. rate: %3.2f;\t%3i %% swap success\t",rank,R,beta,it,acc_rate,swaps);
      mcmc_print_stats(kernel, stdout);
      mcmc_adapt(kernel, acc_rate);
      acc_c = 0;
      swaps=0;
    }
  }
  return EXIT_SUCCESS;
}

int mcmc_foreach(int rank, int R, size_t SampleSize, ode_model_parameters *omp, mcmc_kernel *kernel, hdf5block_t *h5block, void *buffer, main_options *option){
  clock_t ct=clock();
  gsl_matrix *log_para_chunk;
  gsl_vector *log_post_chunk;
  int D=get_number_of_MCMC_variables(omp);
  log_para_chunk=gsl_matrix_alloc(CHUNK,D);
  log_post_chunk=gsl_vector_alloc(CHUNK);
  gsl_vector_view current;
  gsl_vector_view x_state;
  int swaps = 0; // swap success counter
  int acc=no;    // acceptance flag
  int acc_c=0;   // acceptance counter
  double acc_rate=0.0;
  size_t it;
  int master=no;
  int DEST;
  double beta=mcmc_get_beta(kernel);
  herr_t status;
  int resume_EC;
  int last_chunk=no, not_written_yet=yes;
  for (it = 0; it < SampleSize; it++) {
    mcmc_sample(kernel, &acc);
    last_chunk=SampleSize-it<CHUNK;
    if (acc && last_chunk && not_written_yet){
      resume_EC=write_resume_state(option->resume_file, rank, R, kernel);
      if (resume_EC==EXIT_SUCCESS){
	not_written_yet=no;
      }
    }
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
      fprintf(stdout, "# [rank % 2i/% 2i; β=%5f; %3li%% done] (it %5li)\tacc. rate: %3.2f;\t%3i %% swap success\t",rank,R,beta,(100*it)/SampleSize,it,acc_rate,swaps);
      mcmc_print_stats(kernel, stdout);
      acc_c = 0;

      // print chunk to hdf5 file
      status=h5write_current_chunk(h5block,log_para_chunk,log_post_chunk);
      h5block->offset[0]=it+1;
      swaps=0;
    }
  }
  assert(status==0);
  // write remaining data to the output hdf5 file
  int Rest=SampleSize % CHUNK;
  printf("[main] last iteration done %i points remain to write.\n",Rest);
  if (Rest > 0){
    h5block->chunk_size[0]=Rest;
    h5block->chunk_size[1]=D;
    h5block->para_chunk_id=H5Screate_simple(2, h5block->chunk_size, NULL);    
    
    h5block->chunk_size[0]=Rest;
    h5block->chunk_size[1]=1;
    h5block->post_chunk_id=H5Screate_simple(2, h5block->chunk_size, NULL);
    
    printf("[main] writing the remaining %i sampled parametrisations to file.\n",Rest);
    h5block->block[0]=Rest;
    h5block->block[1]=D;
    display_chunk_properties(h5block);
    status = H5Sselect_hyperslab(h5block->para_dataspace_id, H5S_SELECT_SET, h5block->offset, h5block->stride, h5block->count, h5block->block);
    H5Dwrite(h5block->parameter_set_id, H5T_NATIVE_DOUBLE, h5block->para_chunk_id, h5block->para_dataspace_id, H5P_DEFAULT, log_para_chunk->data);
    
    h5block->block[1]=1;
    printf("[main] writing their %i log-posterior values to file.\n",Rest);
    display_chunk_properties(h5block);

    status &= H5Sselect_hyperslab(h5block->post_dataspace_id, H5S_SELECT_SET, h5block->offset, h5block->stride, h5block->count, h5block->block);
    H5Dwrite(h5block->posterior_set_id, H5T_NATIVE_DOUBLE, h5block->post_chunk_id, h5block->post_dataspace_id, H5P_DEFAULT, log_post_chunk->data);
    assert(status>=0);
  }

  // annotate written sample with all necessary information
  //printf("[main] writing some annotation about the sampled points as hdf5 attributes.\n");
  status&=H5LTset_attribute_int(h5block->file_id, "LogParameters", "MPI_RANK", &rank, 1);
  status&=H5LTset_attribute_ulong(h5block->file_id, "LogParameters", "SampleSize", &SampleSize, 1);
  status&=H5LTset_attribute_double(h5block->file_id, "LogParameters", "InverseTemperature_Beta", &beta, 1);
  
  ct=clock()-ct;
  double sampling_time=((double) ct)/((double) CLOCKS_PER_SEC);
  int ts=round(sampling_time);
  int hms[3]; // hours, minutes, seconds
  hms[0]=ts/3600;
  hms[1]=(ts%3600)/60;
  hms[2]=(ts%60);
  printf("# computation time spend sampling: %i:%i:%i\n",hms[0],hms[1],hms[2]);
  
  h5block->size[0]=1;
  status&=H5LTmake_dataset_double (h5block->file_id, "SamplingTime_s", 1, h5block->size, &sampling_time);
  h5block->size[0]=3;
  status&=H5LTmake_dataset_int(h5block->file_id, "SamplingTime_hms", 1, h5block->size, hms);
  
  if(status){
    printf("[rank %i] statistics written to file.\n",rank);
  }
  return EXIT_SUCCESS;
}

herr_t append_meta_properties(hdf5block_t *h5block, double *seed, size_t *BurnInSampleSize, char *h5file, char *lib_base){
  herr_t status;
  status=H5LTset_attribute_string(h5block->file_id, "LogParameters", "ModelLibrary", lib_base);

  status&=H5LTset_attribute_double(h5block->file_id, "LogParameters", "seed", seed, 1);
  status&=H5LTset_attribute_ulong(h5block->file_id, "LogParameters", "BurnIn", BurnInSampleSize, 1);
  status&=H5LTset_attribute_string(h5block->file_id, "LogParameters", "DataFrom", h5file);  
  return status;
}

void print_experiment_information(int rank, int R, ode_model_parameters *omp, gsl_vector *y0){
  int i;
  int C=omp->size->C;
  int N=omp->size->N;
  int NNE=0; // number of normalising experiments
  int LE=0;
  
  for (i=0;i<C;i++){
    if (omp->E[i]->init_y){
      fprintf(stderr,"[main] warning: y(t0) already initialised.\n");
    }else{
      omp->E[i]->init_y=gsl_vector_alloc(N);
      gsl_vector_memcpy(omp->E[i]->init_y,y0);
    }
    if (omp->E[i]->lflag==0) NNE++;
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
}

void display_test_evaluation_results(mcmc_kernel *kernel){
  int i;
  assert(kernel);
  int D=MCMC_DIM(kernel);
  double *x=MCMC_STATE(kernel);
  double *log_p=MCMC_POSTERIOR(kernel);
  printf("# [main] test evaluation of Posterior function done:\n");
  printf("# \tθ=θ₀; LogPosterior(θ|D)=%+g;\n# where θ₀:",log_p[0]); 
  for (i=0;i<D;i++) printf(" %+g ",x[i]);
  printf("\n");
  printf("# [main] LogLikelihood(D|θ):");
  printf("%+g\tLogPrior(θ)=%+g.\n",log_p[1],log_p[2]);
}

main_options get_default_options(char *global_sample_filename_stem, char *lib_name){
  main_options option;
  option.initial_stepsize_rank_factor=1.0;
  option.output_file=global_sample_filename_stem;
  option.library_file=lib_name;
  option.target_acceptance=-0.25;
  option.initial_stepsize =-0.1;
  option.sample_size=-100;
  option.abs_tol=ODE_SOLVER_ABS_ERR;
  option.rel_tol=ODE_SOLVER_REL_ERR;
  return option;
}

int main(int argc, char* argv[]) {
  int i=0;
  int warm_up=0; // sets the number of burn in points at command line
  char lib_name[BUFSZ];
  ode_model_parameters omp[1];
  omp->size=(problem_size*) malloc(sizeof(problem_size));

  char global_sample_filename_stem[BUFSZ]="Sample.h5"; // filename basis
  char rank_sample_file[BUFSZ]; // filename for sample output
  char resume_filename[BUFSZ]="resume.h5";
  double seed = 1;
  double gamma= 2;
  double t0=-1;
  int sampling_action=SMPL_FRESH;
  
  int start_from_prior=no;
  int sensitivity_approximation=no;

  main_options cnf_options=get_default_options(global_sample_filename_stem, lib_name);
  
  MPI_Init(&argc,&argv);
  int rank,R;
  MPI_Comm_size(MPI_COMM_WORLD,&R);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  char *h5file=NULL;

  gsl_set_error_handler_off();
  
  /* process command line arguments
   */  
  for (i=0;i<argc;i++){
    if (strcmp(argv[i],"-p")==0 || strcmp(argv[i],"--prior-start")==0) {
      start_from_prior=1;
    } else if (strcmp(argv[i],"-d")==0 || strcmp(argv[i],"--hdf5")==0) {
      h5file=argv[i+1];
    } else if (strcmp(argv[i],"-t")==0 || strcmp(argv[i],"--init-at-t")==0) {
      t0=strtod(argv[i+1],NULL);
      //printf("[main] t0=%f\n",t0);
    } else if (strcmp(argv[i],"-w")==0 || strcmp(argv[i],"--warm-up")==0) warm_up=strtol(argv[i+1],NULL,10);
    else if (strcmp(argv[i],"--resume")==0 || strcmp(argv[i],"-r")==0) sampling_action=SMPL_RESUME;
    else if (strcmp(argv[i],"--sens-approx")==0) sensitivity_approximation=1;
    else if (strcmp(argv[i],"-l")==0) strcpy(cnf_options.library_file,argv[i+1]);
    //    else if (strcmp(argv[i],"-n")==0) Tuning=0;
    else if (strcmp(argv[i],"-s")==0) cnf_options.sample_size=strtol(argv[i+1],NULL,0);
    else if (strcmp(argv[i],"-o")==0) strncpy(cnf_options.output_file,argv[i+1],BUFSZ);
    else if (strcmp(argv[i],"-a")==0) cnf_options.target_acceptance=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-i")==0 || strcmp(argv[i],"--initial-step-size")==0) cnf_options.initial_stepsize=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-m")==0 || strcmp(argv[i],"--initial-step-size-rank-multiplier")==0) cnf_options.initial_stepsize_rank_factor=strtod(argv[i+1],NULL);

    else if (strcmp(argv[i],"-g")==0) gamma=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"--abs-tol")==0) cnf_options.abs_tol=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"--rel-tol")==0) cnf_options.rel_tol=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"--seed")==0) seed=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--help")==0) {
      print_help();
      MPI_Abort(MPI_COMM_WORLD,0);
    }
  }
  
  seed=seed*137+13*rank;

  /* load Data from hdf5 file
   */
  if (h5file){
    printf("# [main] (rank %i) reading hdf5 file, loading data.\n",rank);
    fflush(stdout);
    read_data(h5file,omp);
    fflush(stdout);
  } else {
    fprintf(stderr,"# [main] (rank %i) no data provided (-d option), exiting.\n",rank);
    MPI_Abort(MPI_COMM_WORLD,-1);
  }

  
  /* load model from shared library
   */
  ode_model *odeModel = ode_model_loadFromFile(lib_name);  /* alloc */
  if (!odeModel) {
    fprintf(stderr, "# [main] (rank %i) Library %s could not be loaded.\n",rank,lib_name);
    exit(1);
  } else printf( "# [main] (rank %i) Library %s loaded.\n",rank, lib_name);
  
  /* construct an output file from rank, library name, and user
   * supplied string.
   */
  char *dot;
  char *lib_base;
  lib_base=basename(lib_name);
  dot=strchr(lib_base,'.');
  dot[0]='\0';
  sprintf(resume_filename,"%s_resume_%02i.h5",lib_base,rank);
  sprintf(rank_sample_file,"mcmc_rank_%02i_of_%i_%s_%s",rank,R,lib_base,basename(cnf_options.output_file));
  cnf_options.output_file=rank_sample_file;
  cnf_options.resume_file=resume_filename;
  
  /* allocate a solver for each experiment for possible parallelization
   */
  ode_solver **solver;
  int c,C=omp->size->C;
  int c_success=0;
  solver=malloc(sizeof(ode_solver*)*C);
  for (c=0;c<C;c++){
    solver[c]=ode_solver_alloc(odeModel);
    if (solver[c]) c_success++;
  }
  if (c_success==C) {
    printf("# [main] Solver[0:%i] for «%s» created.\n",C,lib_base);
  } else {
    fprintf(stderr, "# [main] Solvers for «%s» could not be created.\n",lib_base);
    ode_model_free(odeModel);
    MPI_Abort(MPI_COMM_WORLD,-1);
  }

  /* sensitivity analysis is not feasible for large models. So, it can
   *  be turned off.
   */
  if (sensitivity_approximation){
    //printf("# [main] experimental: Sensitivity approximation activated.\n");
    for (c=0;c<C;c++) ode_solver_disable_sens(solver[c]);
    /* also: make sensitivity function unavailable; that way
     * ode_model_has_sens(model) will return «FALSE»;
     */
    odeModel->vf_sens=NULL;
  }
  
  /* init solver 
   */
  realtype solver_param[3] = {cnf_options.abs_tol, cnf_options.rel_tol, 0};

  const char **x_name=ode_model_get_var_names(odeModel);
  const char **p_name=ode_model_get_param_names(odeModel);
  const char **f_name=ode_model_get_func_names(odeModel);
  
  /* define local variables for parameters and inital conditions */
  int N = ode_model_getN(odeModel);
  int P = ode_model_getP(odeModel);
  int F = ode_model_getF(odeModel);
  
  omp->size->N=N;
  omp->size->P=P;
  omp->size->F=F;
  omp->t0=t0;
  ode_model_parameters_alloc(omp);
  ode_model_parameters_link(omp);
  fflush(stdout);

  /* get default parameters from the model file
   */
  double p[P];
  gsl_vector_view p_view=gsl_vector_view_array(p,P);
  ode_model_get_default_params(odeModel, p, P);
  if (rank==0)  gsl_printf("default parameters",&(p_view.vector),GSL_IS_DOUBLE | GSL_IS_VECTOR);
  omp->solver=solver;

  /* All MCMC meta-parameters (like stepsize) here are positive (to
   * make sense). Some command line arguments can override parameters
   * read from files; but, input files are processed after the command
   * line parameters. So, to check whether default parameters were
   * altered by the command line, the variable declaration defaults
   * are negative at first. Alterations to some meta-parameter p can
   * be checked by: if (cnf_options.p<0)
   * cnf_options.p=read_from_file(SOME FILE);
   */
  cnf_options.initial_stepsize=fabs(cnf_options.initial_stepsize);
  cnf_options.target_acceptance=fabs(cnf_options.target_acceptance);
  cnf_options.sample_size=fabs(cnf_options.sample_size);

  /* load default initial conditions
   */
  double y[N];
  gsl_vector_view y_view=gsl_vector_view_array(y,N);
  ode_model_get_initial_conditions(odeModel, y, N);
  
  print_experiment_information(rank,R,omp,&(y_view.vector));

  /* initialize the ODE solver with initial time t, default ODE
   * parameters p and default initial conditions of the state y; In
   * addition error tolerances are set and sensitivity initialized.
   */
  //printf("# [main] (rank %i) init ivp: t0=%g\n",rank,omp->t0);
  for (c=0;c<C;c++){
    ode_solver_init(solver[c], omp->t0, omp->E[c]->init_y->data, N, p, P);
    //printf("# [main] solver initialised.\n");    
    ode_solver_setErrTol(solver[c], solver_param[1], &solver_param[0], 1);
    if (ode_model_has_sens(odeModel)) {
      ode_solver_init_sens(solver[c], omp->E[0]->yS0->data, P, N);
    }
  }
  /* An smmala_model is a struct that contains the posterior
   * probablity density function and a pointer to its parameters and
   * pre-allocated work-memory.
   */
  smmala_model* model = smmala_model_alloc(LogPosterior, NULL, omp);
  if (model){
    printf("[main] (rank %i) smmala_model allocated.\n",rank);
  }else{
    fprintf(stderr,"[main] (rank %i) smmala_model could not be allocated.\n",rank);
    MPI_Abort(MPI_COMM_WORLD,-1);
  }
  
  /* initial parameter values; after allocating an mcmc_kernel of the
   * right dimensions we set the initial Markov chain state from
   * either the model's default parametrization p, the prior's μ, or the
   * state of a previously completed mcmc run (resume).
   */
  int D=omp->size->D;
  double init_x[D];
  double beta=assign_beta(rank,R,round(gamma));
  double tgac=cnf_options.target_acceptance;
  double m=cnf_options.initial_stepsize_rank_factor;
  double step=cnf_options.initial_stepsize;
  if (m>1.0 && rank>0) step*=gsl_pow_int(m,rank);
  mcmc_kernel* kernel = smmala_kernel_alloc(beta,D,step,model,seed,tgac);
  
  int resume_load_status;
  if (sampling_action==SMPL_RESUME){
    resume_load_status=load_resume_state(resume_filename, rank, R, kernel);
    assert(resume_load_status==EXIT_SUCCESS);
    for (i=0;i<D;i++) init_x[i]=kernel->x[i];
  } else if (start_from_prior){     
    if (rank==0) printf("# [main] setting initial mcmc vector to prior mean.\n");
    for (i=0;i<D;i++) init_x[i]=gsl_vector_get(omp->prior->mu,i);
  } else {
    if (rank==0) printf("# [main] setting mcmc initial value to log(default parameters)\n");
    for (i=0;i<D;i++) init_x[i]=gsl_sf_log(p[i]);
  }
  fflush(stdout);
  //display_prior_information(omp->prior);
  
  /* here we initialize the mcmc_kernel; this makes one test
   * evaluation of the log-posterior density function. 
   */
  
  /*
  if (rank==0){
    printf("# [main] initializing MCMC.\n");
    printf("# [main] init_x:");
    for (i=0;i<D;i++) printf(" %g ",init_x[i]);
    printf("\n");
  }
  */
  
  mcmc_init(kernel, init_x);
  /* display the results of that test evaluation
   *
   */
  if (rank==0){
    printf("# [main] rank %i init complete .\n",rank);
    display_test_evaluation_results(kernel);
    ode_solver_print_stats(solver[0], stdout);
    fflush(stdout);
    fflush(stderr);  
  }
  
  size_t SampleSize = cnf_options.sample_size;  
  
  /* in parallel tempering th echains can swap their positions;
   * this buffers the communication between chains.
   */
  void *buffer=(void *) smmala_comm_buffer_alloc(D);
  
  /* Initialization of burin in length
   */
  size_t BurnInSampleSize;
  if (warm_up==0){
    BurnInSampleSize = 7 * (int) sqrt(cnf_options.sample_size);
  } else {
    BurnInSampleSize=warm_up;
  }
  if (rank==0){
    printf("# Performing Burn-In with step-size (%g) tuning: %lu iterations\n",get_step_size(kernel),BurnInSampleSize);
    fflush(stdout);
  }
  /* Burn In: these iterations are not recorded, but are used to find
   * an acceptable step size for each temperature regime.
   */
  int mcmc_error;
  mcmc_error=burn_in_foreach(rank,R, BurnInSampleSize, omp, kernel, buffer);
  assert(mcmc_error==EXIT_SUCCESS);
  if (rank==0){
    fprintf(stdout, "\n# Burn-in complete, sampling from the posterior.\n");
  }
  /* this struct contains all necessary id's and size arrays
   * for writing sample data to an hdf5 file in chunks
   */
  hdf5block_t *h5block = h5block_init(cnf_options.output_file,
				      omp,SampleSize,
				      x_name,p_name,f_name);
  
  /* The main loop of MCMC sampling
   * these iterations are recorded and saved to an hdf5 file
   * the file is set up and identified via the h5block variable.
   */  
  mcmc_error=mcmc_foreach(rank, R, SampleSize, omp, kernel, h5block, buffer, &cnf_options);
  assert(mcmc_error==EXIT_SUCCESS);
  append_meta_properties(h5block,&seed,&BurnInSampleSize, h5file, lib_base);
  h5block_close(h5block);

  /* clear memory */
  smmala_model_free(model);
  mcmc_free(kernel);
  ode_model_parameters_free(omp);
  MPI_Finalize();
  return EXIT_SUCCESS;
}
