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
#include "normalisation_sd.h"
#include "../mcmc/smmala.h"
#include "../ode/ode_model.h"
#include "../mcmc/smmala_posterior.h"
#include "diagnosis_output.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "omp.h"
#include "options.h"
#include "../mcmc/model_parameters_smmala.h"
// define target block types
#define INTEGER_BLOCK 1
#define DOUBLE_BLOCK 2

#define yes 1
#define no 0
#define CHUNK 100
// sampling actions 
#define SMPL_RESUME 1
#define SMPL_FRESH 0
#define SMPL_RESUME_TUNE 2

#define BUFSZ 2048


/* collects all parameters and size arrays needed for hdf5 functions
 */
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

/* closes all hdf5 id's and frees h5block
 */
int /*always returns success*/ h5block_close(hdf5block_t *h5block){
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

/* initializes hdf5 struct and writes some of the initially known
   problem properties such as state variable names into hdf5 file*/
hdf5block_t* /*freshly allocated struct with ids and size parameters*/
h5block_init(char *output_file, /*will create this file for writing*/
	     ode_model_parameters *omp, /*contains MCMC problem description*/
	     size_t Samples, /* MCMC sample size*/
	     const char **x_name, /*ODE model State Variable names, array of strings*/
	     const char **p_name, /*ODE model parameter names, array of strings*/
	     const char **f_name)/*ODE model output function names, array of strings*/
{
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
  NameWriteError|=H5LTmake_dataset_string(file_id,"StateVariableNames",x_names);
  NameWriteError|=H5LTmake_dataset_string(file_id,"ParameterNames",p_names);
  NameWriteError|=H5LTmake_dataset_string(file_id,"OutputFunctionNames",f_names);
  if (NameWriteError){
    fprintf(stderr,"[h5block_init] writing (x,p,f)-names into hdf5 file failed.");
  }
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

/*assigns a temperature @code beta@ to an MPI rank*/
double /*beta*/ assign_beta(int rank,/*MPI rank*/ int R, /*MPI Comm size*/ int gamma)/*exponent: @code beta=(1-rank/R)^gamma@*/{
  double x=(double)(R-rank)/(double) R;
  double b=-1;
  assert(gamma>=1);
  b=gsl_pow_int(x, gamma);
  assert(b>=0 && b<=1.0);
  return b;
}

/*output of @code ./ode_smmala --help@*/
void print_help(){
  printf("Usage ($SOMETHING are values you choose, written as bash variables):\n");
  printf("-a $ACCEPTANCE_RATE\n");
  printf("\t\t\tTarget acceptance value (all markov chains will be tuned for this acceptance).\n\n");
  printf("-d, --hdf5 ./data.h5\n");
  printf("\t\t\tdata.h5 is a file that contains the data points and the conditions of measurement in hdf5 format. A suitable h5 file is produced by the hdf5_import program bundled with ode_smmala.\n\n");
  printf("-g $G\n");
  printf("\t\t\tThis will define how the inverse MCMC temperatures β are chosen: β = (1-rank/R)^G, where R is MPI_Comm_Size.\n\n");
  printf("-i $STEP_SIZE\n");
  printf("\t\t\tThe initial step size of each markov chain, this will usually be tuned later to get the desired acceptance rate $A (-a $A).\n\n");
  printf("-l ./ode_model.so\n");
  printf("\t\t\tode_model.so is a shared library containing the CVODE functions of the model.\n\n");
  printf("-m $M\n");
  printf("\t\t\tIf this number is larger than 1.0, each MPI rank will get a different initial step size s: step_size(rank)=STEP_SIZE*M^(rank).\n\n");
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

/* pseudo graphical display of sample chunk statistical properties,
   only useful with one MPI worker as it prints to terminal.
 */
void print_chunk_graph(gsl_matrix *X,/*sub-sample of CHUNK rows*/ gsl_vector *lP)/*log-Posterior values, unnormalized*/{
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

/*debug function: shows whether the size, offset, count and stride have been set to sensible values*/
void display_chunk_properties(hdf5block_t *h5block)/*structure holding the hdf5 parameters*/{
  hsize_t *offset=h5block->offset;
  hsize_t *block=h5block->block;
  hsize_t *count=h5block->count;
  hsize_t *stride=h5block->stride;
  printf("[%s]",__func__);
  printf("  offset: %lli×%lli;\n",offset[0],offset[1]);
  printf("   block: %lli×%lli;\n",block[0],block[1]);
  printf("  stride: %lli×%lli;\n",stride[0],stride[1]);
  printf("   count: %lli×%lli.\n",count[0],count[1]);
}

/*writes a sampled chunk into the appropriate hyperslab of hdf5 file*/
herr_t /*hdf5 error type*/
h5write_current_chunk
(hdf5block_t *h5block,/*holds hdf5 properties and ids*/
 gsl_matrix *log_para_chunk, /*log-parameter chunk*/
 gsl_vector *log_post_chunk)/*log-posterior value chunk*/{
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
/*this executes a sampling loop, without recording the values, it adapts step size to get target acceptance*/
int /*always returns success*/
burn_in_foreach
(int rank, /*MPI rank*/
 int R, /*MPI Comm size*/
 size_t BurnInSampleSize, /*number of iterations for tuning*/
 ode_model_parameters *omp, /*ODE problem definition, allocated space*/
 mcmc_kernel *kernel, /*MCMC kernel struct*/
 void *buffer)/*MPI communication buffer, a deep copy of @code kernel@*/{
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
      fflush(stdout);
      mcmc_adapt(kernel, acc_rate);
      acc_c = 0;
      swaps=0;
    }
  }
  return EXIT_SUCCESS;
}

/*main mcmc loop, records sampled values, performs no tuning of the step size.*/
int /*always returns success*/
mcmc_foreach
(int rank, /*MPI rank*/
 int R, /*MPI Comm size*/
 size_t SampleSize, /*number of iterations of MCMC*/
 ode_model_parameters *omp, /*ODE problem cpecification and pre-allocated space*/
 mcmc_kernel *kernel, /* MCMC kernel sturct, holds MCMC-algorithm's parameters*/
 hdf5block_t *h5block, /*defines the hdf5 file to write into, holds ids and sizes*/
 void *buffer, /* for MPI communication, similar to kernel*/
 main_options *option)/*options from defaults, files and command line*/{
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
      fflush(stdout);
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
  printf("[%s] last iteration done %i points remain to write.\n",__func__,Rest);
  if (Rest > 0){
    h5block->chunk_size[0]=Rest;
    h5block->chunk_size[1]=D;
    h5block->para_chunk_id=H5Screate_simple(2, h5block->chunk_size, NULL);    
    
    h5block->chunk_size[0]=Rest;
    h5block->chunk_size[1]=1;
    h5block->post_chunk_id=H5Screate_simple(2, h5block->chunk_size, NULL);
    
    printf("[%s] writing the remaining %i sampled parametrisations to file.\n",__func__,Rest);
    h5block->block[0]=Rest;
    h5block->block[1]=D;
    display_chunk_properties(h5block);
    status = H5Sselect_hyperslab(h5block->para_dataspace_id, H5S_SELECT_SET, h5block->offset, h5block->stride, h5block->count, h5block->block);
    H5Dwrite(h5block->parameter_set_id, H5T_NATIVE_DOUBLE, h5block->para_chunk_id, h5block->para_dataspace_id, H5P_DEFAULT, log_para_chunk->data);
    
    h5block->block[1]=1;
    printf("[%s] writing their %i log-posterior values to file.\n",__func__,Rest);
    display_chunk_properties(h5block);

    status |= H5Sselect_hyperslab(h5block->post_dataspace_id, H5S_SELECT_SET, h5block->offset, h5block->stride, h5block->count, h5block->block);
    H5Dwrite(h5block->posterior_set_id, H5T_NATIVE_DOUBLE, h5block->post_chunk_id, h5block->post_dataspace_id, H5P_DEFAULT, log_post_chunk->data);
    assert(status>=0);
  }

  // annotate written sample with all necessary information
  //printf("[main] writing some annotation about the sampled points as hdf5 attributes.\n");
  status|=H5LTset_attribute_int(h5block->file_id, "LogParameters", "MPI_RANK", &rank, 1);
  status|=H5LTset_attribute_ulong(h5block->file_id, "LogParameters", "SampleSize", &SampleSize, 1);
  status|=H5LTset_attribute_double(h5block->file_id, "LogParameters", "InverseTemperature_Beta", &beta, 1);
  
  ct=clock()-ct;
  double sampling_time=((double) ct)/((double) CLOCKS_PER_SEC);
  int ts=round(sampling_time);
  int hms[3]; // hours, minutes, seconds
  hms[0]=ts/3600;
  hms[1]=(ts%3600)/60;
  hms[2]=(ts%60);
  printf("# computation time spend sampling: %i:%i:%i\n",hms[0],hms[1],hms[2]);
  
  h5block->size[0]=1;
  status|=H5LTmake_dataset_double (h5block->file_id, "SamplingTime_s", 1, h5block->size, &sampling_time);
  h5block->size[0]=3;
  status|=H5LTmake_dataset_int(h5block->file_id, "SamplingTime_hms", 1, h5block->size, hms);
  
  if(status){
    printf("[rank %i] statistics written to file.\n",rank);
  }
  return EXIT_SUCCESS;
}

/*writes properties of the current run related to implementation and command line choices*/
herr_t /*hdf5 error*/
append_meta_properties(hdf5block_t *h5block,/*hdf5 file ids*/
		       double *seed,/*random number seed*/
		       size_t *BurnInSampleSize, /*tuning iterations*/
		       char *h5file, /*name of hdf5 file containing the experimental data and prior set-up*/
		       char *lib_base)/*basename of the library file @code .so@ file*/{
  herr_t status;
  int omp_n,omp_np,i=1;
  status=H5LTset_attribute_double(h5block->file_id, "LogParameters", "seed", seed, 1);
  status|=H5LTset_attribute_ulong(h5block->file_id, "LogParameters", "BurnIn", BurnInSampleSize, 1);
  status|=H5LTset_attribute_string(h5block->file_id, "LogParameters", "DataFrom", h5file);
  status|=H5LTmake_dataset_string(h5block->file_id,"Model",lib_base);
  // here we make a short test to see what the automatic choice of the
  // number of threads turns out to be.
#pragma omp parallel private(omp_n,omp_np) reduction(+:i)
  {
    i=1;
    omp_n=omp_get_num_threads();
    omp_np=omp_get_num_procs();
  }
  if (i!=omp_n){
    fprintf(stderr,"[append_meta_properties] warning: finding out number of threads possibly failed reduction of (n×1: %i) != get_num_threads():%i.\n",i,omp_n);
  } 
  h5block->size[0]=1;
  h5block->size[1]=1;
  status|=H5LTmake_dataset_int(h5block->file_id,"OMP_NUM_THREADS",1,h5block->size,&omp_n);
  status|=H5LTmake_dataset_int(h5block->file_id,"OMP_NUM_PROCS",1,h5block->size,&omp_np);
  return status;
}

/*prints how many experiments are normalization experiments and sets
  each experiment's initial conditions if not previously set*/
void print_experiment_information
(int rank,/*MPI rank*/
 int R, /*MPI Comm size*/
 ode_model_parameters *omp, /*ODE model parameters*/
 gsl_vector *y0)/*globally set initial conditions*/{
  int i;
  int C=omp->size->C;
  int N=omp->size->N;
  int NNE=0; // number of normalising experiments
  int LE=0;
  
  for (i=0;i<C;i++){
    if (omp->E[i]->init_y){
      fprintf(stderr,"[%s] y(t0) set by h5 file in Experiment %i of %i. I'll leave it unchanged\n",__func__,i,C);
    }else{
      fprintf(stderr,"[%s] y(t0) is not yet set for Experiment %i of %i. I'll use the global default from the model's «.so» file (set in vfgen file).\n",__func__,i,C);
      omp->E[i]->init_y=gsl_vector_alloc(N);
      gsl_vector_memcpy(omp->E[i]->init_y,y0);
    }
    if (omp->E[i]->lflag==0) NNE++;
  }
  fflush(stdout);
  LE=C-NNE;
  if (rank==0){
    printf("[%s] There are %i experiments",__func__,C);
    if (NNE>0){
      printf(", %i of which ",NNE);
      if (NNE==1) printf("is");
      else printf("are");
      printf(" used only for the normalisation of the %i experiments that explicitly contribute to the LogLikelihood(NormalisedData[1:%i]|θ).\n",LE,LE); 	
    }  else printf(".\n");
  }
  fflush(stdout);
}

/*Kernel init makes a test evakuation of log-posterior pdf; this function prints the results (uses a couple of unicode characters)*/
void display_test_evaluation_results(mcmc_kernel *kernel)/*MCMC kernel struct*/{
  int i;
  assert(kernel);
  int D=MCMC_DIM(kernel);
  double *x=MCMC_STATE(kernel);
  double *log_p=MCMC_POSTERIOR(kernel);
  printf("# [%s] test evaluation of Posterior function done:\n",__func__);
  printf("# \tθ=θ₀; LogPosterior(θ|D)=%+g;\n# where θ₀:",log_p[0]); 
  for (i=0;i<D;i++) printf(" %+g ",x[i]);
  printf("\n");
  printf("[%s] LogLikelihood(D|θ):",__func__);
  printf("%+g\tLogPrior(θ)=%+g.\n",log_p[1],log_p[2]);
}

/* Calculates the normalisation constant of the likelihood function:
 * 1/sqrt(2*pi*sigma²)^beta for all data points (product).
 * This is done in log-space
 */
int /*error flag*/
pdf_normalisation_constant(ode_model_parameters *omp)/*pre-allocated storage for simulation results, used in LogLikelihood calculations*/{
  assert(omp && omp->size);
  int c,C=get_number_of_experimental_conditions(omp);
  int i,F=get_number_of_model_outputs(omp);
  int t,T;
  double E_lN,lN=0; // log normalisation constant;
  double stdv;
  for (c=0;c<C;c++){
    T=omp->E[c]->t->size;
    E_lN=-0.5*(M_LN2+M_LNPI);
    for (t=0;t<T;t++){
      for (i=0;i<F;i++){
	stdv=gsl_vector_get(omp->E[c]->sd_data[t],i);
	if (gsl_finite(stdv)){
	  E_lN-=gsl_sf_log(stdv);
	}
      }
    }
    omp->E[c]->pdf_lognorm=E_lN;
    lN+=E_lN;
  }
  omp->pdf_lognorm=lN;
  return EXIT_SUCCESS;
}


/* Initializes MPI,
 * loads defaults, 
 *       command line arguments,
 *       hdf5 data,
 *       ode model from shared library @code dlopen@
 * allocates kernel, 
 *           ode model parameters
 *           MPI communivcation buffers
 * calls MCMC routines
 * finalizes and frees (most) structs
 */
int/*always returns success*/
main(int argc,/*count*/ char* argv[])/*array of strings*/ {
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

  main_options *cnf_options=default_options(global_sample_filename_stem, lib_name);
  
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
    } else if (strcmp(argv[i],"-w")==0 || strcmp(argv[i],"--warm-up")==0) warm_up=strtol(argv[i+1],NULL,10);
    else if (strcmp(argv[i],"--resume")==0 || strcmp(argv[i],"-r")==0) sampling_action=SMPL_RESUME;
    else if (strcmp(argv[i],"--sens-approx")==0) sensitivity_approximation=1;
    else if (strcmp(argv[i],"-l")==0) strcpy(cnf_options->library_file,argv[i+1]);
    //    else if (strcmp(argv[i],"-n")==0) Tuning=0;
    else if (strcmp(argv[i],"-s")==0) cnf_options->sample_size=strtol(argv[i+1],NULL,0);
    else if (strcmp(argv[i],"-o")==0) strncpy(cnf_options->output_file,argv[i+1],BUFSZ);
    else if (strcmp(argv[i],"-a")==0) cnf_options->target_acceptance=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-i")==0 || strcmp(argv[i],"--initial-step-size")==0) cnf_options->initial_stepsize=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-m")==0 || strcmp(argv[i],"--initial-step-size-rank-multiplier")==0) cnf_options->initial_stepsize_rank_factor=strtod(argv[i+1],NULL);

    else if (strcmp(argv[i],"-g")==0) gamma=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"--abs-tol")==0) read_abs_tol(cnf_options->abs_tol,argv[i+1]);
    else if (strcmp(argv[i],"--rel-tol")==0) cnf_options->rel_tol=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"--seed")==0) seed=strtod(argv[i+1],NULL);
    else if (strcmp(argv[i],"-h")==0 || strcmp(argv[i],"--help")==0) {
      print_help();
      MPI_Abort(MPI_COMM_WORLD,0);
    }
  }
  
  seed=seed*137+13*rank;

  /* load model from shared library
   */
  ode_model *odeModel = ode_model_loadFromFile(lib_name);  /* alloc */
  if (!odeModel) {
    fprintf(stderr, "[%s] (rank %i) Library %s could not be loaded.\n",__func__,rank,lib_name);
    exit(1);
  } else printf( "[%s] (rank %i) Library %s loaded.\n",__func__,rank, lib_name);
  
  /* construct an output file from rank, library name, and user
   * supplied string.
   */
  char *dot;
  char *lib_base;
  lib_base=basename(lib_name);
  dot=strchr(lib_base,'.');
  if (dot)  dot[0]='\0';
  sprintf(resume_filename,"%s_resume_%02i.h5",lib_base,rank);
  sprintf(rank_sample_file,"mcmc_rank_%02i_of_%i_%s_%s",rank,R,lib_base,basename(cnf_options->output_file));
  cnf_options->output_file=rank_sample_file;
  cnf_options->resume_file=resume_filename;

  /* local variables for parameters and inital conditions as presented
     in ode model lib: */
  int N = ode_model_getN(odeModel);
  int P = ode_model_getP(odeModel);
  int F = ode_model_getF(odeModel);
  printf("[%s] problem size: %i state variables, %i parameters [coefficients and inputs], %i output functions.\n",__func__,N,P,F);
  const char **x_name=ode_model_get_var_names(odeModel);
  const char **p_name=ode_model_get_param_names(odeModel);
  const char **f_name=ode_model_get_func_names(odeModel);
  omp->model_x.name=x_name;
  omp->model_x.size=N;
  omp->model_p.name=p_name;
  omp->model_p.size=P;
  omp->model_f.name=f_name;
  omp->model_f.size=F;

  /* load Data from hdf5 file
   */
  if (h5file){
    printf("[%s] (rank %i) reading hdf5 file, loading data.\n",__func__,rank);
    fflush(stdout);
    read_data(h5file,omp);
    fflush(stdout);
  } else {
    fprintf(stderr,"[%s] (rank %i) no data provided (-d option), exiting.\n",__func__,rank);
    MPI_Abort(MPI_COMM_WORLD,-1);
  }
  fflush(stdout);
  /* allocate a solver for each experiment for parallelization
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
    printf("[%s] Solver[0:%i] for «%s» created.\n",__func__,C,lib_base);
  } else {
    fprintf(stderr, "[%s] Solvers for «%s» could not be created.\n",__func__,lib_base);
    ode_model_free(odeModel);
    MPI_Abort(MPI_COMM_WORLD,-1);
  }
  fflush(stdout);
  /* sensitivity analysis is not feasible for large models. So, it can
   *  be turned off.
   */
  if (sensitivity_approximation){
    for (c=0;c<C;c++) ode_solver_disable_sens(solver[c]);
    /* also: make sensitivity function unavailable; that way
     * ode_model_has_sens(model) will return «FALSE»;
     */
    odeModel->vf_sens=NULL;
    printf("[%s] sensitivity approximation is turned on.\n",__func__);
  }
  fflush(stdout);
  omp->t0=t0;
    /* save in ode model parameter struct: */
  set_number_of_state_variables(omp,N);
  set_number_of_model_parameters(omp,P);
  set_number_of_model_outputs(omp,F);

  /* ode model parameter struct has pointers for sim results that need
     memory allocation: */
  ode_model_parameters_alloc(omp);
  /* data matrix row views are made */
  ode_model_parameters_link(omp);
  /* necessity of normalisation will be checked and noted for later: */
  data_normalisation(omp);
  //printf("[%s] data normalization done.\n",__func__);
  //for (c=0;c<C;c++){
  //  printf("[%s] E%i needs normalization: %i\n",__func__,c,NEEDS_NORMALISATION(omp->E[c]));
  //}
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
   * be checked by: if (cnf_options->p<0)
   * cnf_options->p=read_from_file(SOME FILE);
   */
  cnf_options->initial_stepsize=fabs(cnf_options->initial_stepsize);
  cnf_options->target_acceptance=fabs(cnf_options->target_acceptance);
  cnf_options->sample_size=fabs(cnf_options->sample_size);

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
  for (c=0;c<C;c++){
    ode_solver_init(solver[c], omp->t0, omp->E[c]->init_y->data, N, p, P);
    ode_solver_setErrTol
      (solver[c],
       cnf_options->rel_tol,
       cnf_options->abs_tol->value, cnf_options->abs_tol->size);
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
    printf("[%s] (rank %i) smmala_model allocated.\n",__func__,rank);
  }else{
    fprintf(stderr,"[%s] (rank %i) smmala_model could not be allocated.\n",__func__,rank);
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
  double tgac=cnf_options->target_acceptance;
  double m=cnf_options->initial_stepsize_rank_factor;
  double step=cnf_options->initial_stepsize;
  if (m>1.0 && rank>0) step*=gsl_pow_int(m,rank);
  pdf_normalisation_constant(omp);
  printf("[%s] (rank %i) likelihood log(normalisation constant): %g\n",__func__,rank,omp->pdf_lognorm);
  mcmc_kernel* kernel = smmala_kernel_alloc(beta,D,step,model,seed,tgac);
  
  int resume_load_status;
  if (sampling_action==SMPL_RESUME){
    resume_load_status=load_resume_state(resume_filename, rank, R, kernel);
    assert(resume_load_status==EXIT_SUCCESS);
    for (i=0;i<D;i++) init_x[i]=kernel->x[i];
  } else if (start_from_prior){     
    if (rank==0) printf("[%s] setting initial mcmc vector to prior mean (size %i).\n",__func__,D);
    for (i=0;i<D;i++) init_x[i]=gsl_vector_get(omp->prior->mu,i);
  } else {
    if (rank==0) printf("[%s] setting mcmc initial value to log(default parameters)\n",__func__);
    for (i=0;i<D;i++) init_x[i]=gsl_sf_log(p[i]);
  }
  
  //display_prior_information(omp->prior);
  fflush(stdout);
  /* here we initialize the mcmc_kernel; this makes one test
   * evaluation of the log-posterior density function. 
   */
  mcmc_init(kernel, init_x);
  /* display the results of that test evaluation
   *
   */
  fflush(stdout);
  if (rank==0){
    printf("[%s] rank %i init complete .\n",__func__,rank);
    fflush(stdout);
    display_test_evaluation_results(kernel);
    ode_solver_print_stats(solver[0], stdout);
    fflush(stdout);
    fflush(stderr);  
  }
  
  size_t SampleSize = cnf_options->sample_size;  
  
  /* in parallel tempering th echains can swap their positions;
   * this buffers the communication between chains.
   */
  void *buffer=(void *) smmala_comm_buffer_alloc(D);
  
  /* Initialization of burin in length
   */
  size_t BurnInSampleSize;
  if (warm_up==0){
    BurnInSampleSize = 7 * (int) sqrt(cnf_options->sample_size);
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
  hdf5block_t *h5block = h5block_init(cnf_options->output_file,
				      omp,SampleSize,
				      x_name,p_name,f_name);
  
  /* The main loop of MCMC sampling
   * these iterations are recorded and saved to an hdf5 file
   * the file is set up and identified via the h5block variable.
   */  
  mcmc_error=mcmc_foreach(rank, R, SampleSize, omp, kernel, h5block, buffer, cnf_options);
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
