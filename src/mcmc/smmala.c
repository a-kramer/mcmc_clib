/*
 * SMMALA implementation
 */
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_linalg.h>
#include <string.h>
#include <mpi.h>
#include "smmala.h"
#include "mv_norm.h"
#include "model_parameters_smmala.h"
#include "../app/diagnosis_output.h"
#include "hdf5.h"
#include "hdf5_hl.h"


#define SWAP(a, b, tmp)  tmp = a; a = b; b = tmp

// fx = beta*lx + px
typedef struct{
  double beta;
  gsl_vector *x;
  double *fx; // LogPosterior(x) = beta * lx + px; fx[0]=beta*fx[1] + fx[2]
  double stepsize;	/* store current step size */
  gsl_vector **dfx;	/* store current log-posterior gradient wrt to x: dfx[i]=beta*dlx[i] + dpx (likelihood term and prior term)*/
  gsl_matrix **Hfx; // current cholesky factor of Fisher Information (FI): beta^2 * Hlx_l + Hpx
  gsl_vector *new_x;		/* tmp working space for new x */
  gsl_vector **new_dfx;
  gsl_matrix **new_Hfx;
  gsl_vector *mean_vec;
  gsl_matrix *cholH_mat;
  double target_acceptance;
} smmala_params;

static smmala_params* smmala_params_alloc(const double beta, const int N, double step_size, const double target_acceptance);

void* smmala_comm_buffer_alloc(int D){
  smmala_params *smmala_buffer=smmala_params_alloc(0,D,0,0);;
  if (smmala_buffer==NULL) {
    perror("could not create communication buffer.\n");
    exit(-1);
  }
  return smmala_buffer;
}

double get_step_size(const mcmc_kernel *kernel){
  smmala_params *state = kernel->kernel_params;
  return state->stepsize;
}

int smmala_exchange_information(mcmc_kernel* kernel, const int DEST, void *buffer){
  MPI_Status status; // MPI_Status contains: MPI_SOURCE, MPI_TAG, MPI_ERROR
  int N=kernel->N; // size of x, dfx; N*N is size of FI;
  smmala_params *state = (smmala_params*) kernel->kernel_params;
  smmala_params *scbuf = buffer;
  int TAG=0;
  int SRC=DEST; // send to and receive from the same process

  //error handling for nerr communication steps
  int nerr=6;
  int i, ec=MPI_SUCCESS, overall_error=MPI_SUCCESS;
  int len, error_class;
  char error_string[MPI_MAX_ERROR_STRING];

  //things to be communicated:
  double *send_buffer[]={state->x->data, state->Hfx[1]->data, state->Hfx[2]->data, state->fx, state->dfx[1]->data, state->dfx[2]->data};
  double *recv_buffer[]={scbuf->x->data, scbuf->Hfx[1]->data, scbuf->Hfx[2]->data, scbuf->fx, scbuf->dfx[1]->data, scbuf->dfx[2]->data};
  int            size[]={             N,                 N*N,                 N*N,         3,                   N,                   N};
  /* int n=sizeof(send_buffer); */
  /* if (n!=nerr) { */
  /*   fprintf(stderr,"sizeof operator doesn't work as expected on arrays.\n"); */
  /*   MPI_Abort(MPI_COMM_WORLD,-1); */
  /* } */
  //int r;
  //MPI_Comm_rank(MPI_COMM_WORLD,&r);
  //printf("[rank %i] exchange information\n",r);  
  // every kernel needs to know x and fx;
  for (i=0;i<nerr;i++){
    /* if (send_buffer[i]==NULL) { */
    /*   printf("[rank %i] send_buffer[%i] is a NULL pointer\n",r,i); */
    /*   MPI_Abort(MPI_COMM_WORLD, overall_error); */
    /* } */
    /* if (recv_buffer[i]==NULL) { */
    /*   printf("[rank %i] recv_buffer[%i] is a NULL pointer\n",r,i); */
    /*   MPI_Abort(MPI_COMM_WORLD, overall_error); */
    /* } */

    ec=MPI_Sendrecv(   send_buffer[i], size[i], MPI_DOUBLE, DEST, TAG,
		       recv_buffer[i], size[i], MPI_DOUBLE, SRC, TAG,
		       MPI_COMM_WORLD, &status);
    overall_error&=ec;
    overall_error&=status.MPI_ERROR;
    if (ec != MPI_SUCCESS){
      MPI_Error_string(ec, error_string, &len);
      fprintf(stderr, "[Comm with rank %3i, part %i] error: %s\n", DEST, i, error_string);
      
      MPI_Error_class(ec, &error_class);
      MPI_Error_string(error_class, error_string, &len);
      fprintf(stderr, "[Comm with rank %3i, part %i] class: %s\n", DEST, i, error_string);
    }
  }
  //printf("[rank %i] exchange information done\n",r);  
  if (overall_error != MPI_SUCCESS){
   fprintf(stderr, "[Comm with rank %3i] Some Communication error occured.\n", DEST);
   MPI_Error_string(overall_error, error_string, &len);
   fprintf(stderr, "[Comm with rank %3i] %s\n", DEST, error_string);
   MPI_Abort(MPI_COMM_WORLD, overall_error);
  }
  //FILE *f;
  //char name[128];
  /* sprintf(name,"rank_%i_comm_results.txt",r); */
  /* f=fopen(name,"a"); */
  /* MPI_Barrier(MPI_COMM_WORLD); */
  /* fprintf(f,"[rank %i] my beta: %+g",r,state->beta); */
  /* fprintf(f,"[rank %i] my x: ",r);   for (i=0;i<N;i++)   fprintf(f," %+g ",state->x->data[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] my dfx1: ",r); for (i=0;i<N;i++)   fprintf(f," %+g ",state->dfx[1]->data[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] my dfx2: ",r); for (i=0;i<N;i++)   fprintf(f," %+g ",state->dfx[2]->data[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] my Hfx1: ",r); for (i=0;i<N*N;i++) fprintf(f," %+g ",state->Hfx[1]->data[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] my Hfx2: ",r); for (i=0;i<N*N;i++) fprintf(f," %+g ",state->Hfx[2]->data[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] my beta*lx+px: ",r);              fprintf(f,"beta × %+g + %g",state->fx[1],state->fx[2]); fprintf(f,"\n"); */
  
  /* MPI_Barrier(MPI_COMM_WORLD); */
  /* fprintf(f,"[rank %i] recv x: ",r);   for (i=0;i<N;i++)   fprintf(f," %+g ",scbuf->x->data[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] recv dfx1: ",r); for (i=0;i<N;i++)   fprintf(f," %+g ",scbuf->dfx[1]->data[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] recv dfx2: ",r); for (i=0;i<N;i++)   fprintf(f," %+g ",scbuf->dfx[2]->data[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] recv Hfx1: ",r); for (i=0;i<N*N;i++) fprintf(f," %+g ",scbuf->Hfx[1]->data[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] recv Hfx2: ",r); for (i=0;i<N*N;i++) fprintf(f," %+g ",scbuf->Hfx[2]->data[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] recv beta*lx+px: ",r);                      fprintf(f,"beta × %+g + %g",scbuf->fx[1],scbuf->fx[2]); fprintf(f,"\n"); */
  /* fclose(f); */
  return overall_error;
}

int smmala_swap_chains(mcmc_kernel* kernel, const int master, const int rank, const int other_rank, void *buffer){
  double a;
  int swap_accepted=0;
  //int N=kernel->N; // size of x, dfx; N*N is size of FI;
  gsl_rng *rng = (gsl_rng*) kernel->rng;
  smmala_params *state = (smmala_params*) kernel->kernel_params;
  smmala_params *recv_buf=buffer;
  int i;
  int TAG=0;
  double beta=state->beta;
  double their_beta;
  MPI_Sendrecv(&beta, 1, MPI_DOUBLE, other_rank, TAG, &their_beta, 1, MPI_DOUBLE, other_rank, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  TAG++;
  int gsl_status;
  gsl_sf_result result;
  double beta_state = (their_beta - beta)*(state->fx[1] - recv_buf->fx[1]);
  
  gsl_status=gsl_sf_exp_e(beta_state, &result);
  if (gsl_status==GSL_SUCCESS && result.err<result.val){
    a=result.val;
  } else {
    //fprintf(stderr,"[smmala_swap_chains] gsl_sf_exp(%f)=%f±%f failed.\n",beta_state,result.val,result.err);
    //fprintf(stderr,"\t(rank %i) swap probability: exp((%+g - %+g)*(%+g - %+g))=%+g\n",rank,their_beta,beta,state->fx[1],recv_buf->fx[1],a);
    a=0.0;
  }
  
  if (master){
    double r1=gsl_rng_uniform(rng);
    swap_accepted=r1<a?1:0;
    //printf("[rank %i (master)] swap accepted: %i\n",rank,swap_accepted);
    MPI_Send(&swap_accepted, 1, MPI_INT, other_rank,TAG,MPI_COMM_WORLD);
  }else{
    MPI_Recv(&swap_accepted, 1, MPI_INT, other_rank,TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (swap_accepted){ // both sides copy their received communication buffer into the kernel
    gsl_vector_memcpy(state->x,recv_buf->x);
    memcpy(state->fx,recv_buf->fx,sizeof(double)*3);
    state->fx[0]=beta*state->fx[1] + state->fx[2];
    for (i=1;i<3;i++) {
      gsl_vector_memcpy(state->dfx[i],recv_buf->dfx[i]);
      gsl_matrix_memcpy(state->Hfx[i],recv_buf->Hfx[i]);
    }
    gsl_vector_memcpy(state->dfx[0],state->dfx[1]);
    gsl_vector_scale(state->dfx[0],beta);
    gsl_vector_add(state->dfx[0],state->dfx[2]);
    //gsl_printf("dfx0",state->dfx[0],GSL_IS_VECTOR|GSL_IS_DOUBLE);

    gsl_matrix_memcpy(state->Hfx[0],state->Hfx[1]);
    gsl_matrix_scale(state->Hfx[0],beta*beta);
    gsl_matrix_add(state->Hfx[0],state->Hfx[2]);
    if (gsl_linalg_cholesky_decomp(state->Hfx[0])!=GSL_SUCCESS) perror("Hfx is not positive definite.");
  }
  return swap_accepted;
}

/* This function writes the state of the markov chain into an hdf5
 * file, which can later be loaded to resume sampling.
 */
int write_resume_state(const char *file_name, int rank, int R, const mcmc_kernel *kernel){
  assert(rank<R);
  assert(file_name);
  assert(kernel);
  hid_t file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hsize_t size[i_prior];
  herr_t h5err=0;
  double *buffer;
  smmala_params *state = kernel->kernel_params;
  assert(state);
  
  // beta; the relaxation parameter of parallel tempering
  size[0]=1;
  buffer=&(state->beta);
  h5err&=H5LTmake_dataset_double(file_id,"beta", 1, size, buffer);
  // step size;
  buffer=&(state->stepsize);
  h5err&=H5LTmake_dataset_double(file_id,"stepsize", 1, size, buffer);
  // MPI rank
  h5err&=H5LTmake_dataset_int(file_id,"MPI_Comm_rank", 1, size, &rank);
  h5err&=H5LTset_attribute_int(file_id,"MPI_Comm_rank", "MPI_Comm_size", &R,1);
  // posterior value fx
  // (likelihood, prior)
  size[0]=3;
  buffer=state->fx;
  h5err&=H5LTmake_dataset_double(file_id,"Posterior", 1, size, buffer);
  h5err&=H5LTset_attribute_string(file_id,"Posterior", "key", "LogPosterior; LogLikelihood; LogPrior");
  h5err&=H5LTset_attribute_string(file_id,"Posterior", "info", "LogPosterior=beta*LogLikelihood+LogPrior");
  assert(h5err==0);
  // state of the Markov chain x
  size[0]=state->x->size;
  buffer=gsl_vector_ptr(state->x,0);
  h5err&=H5LTmake_dataset_double(file_id,"MarkovChainState", 1, size, buffer);
  // Gradients 
  size[0]=state->dfx[i_posterior]->size;
  buffer=gsl_vector_ptr(state->dfx[i_posterior],0);
  h5err&=H5LTmake_dataset_double(file_id,"LogPosteriorGradient", 1, size, buffer);
  buffer=gsl_vector_ptr(state->dfx[i_likelihood],0);
  h5err&=H5LTmake_dataset_double(file_id,"LogLikelihoodGradient", 1, size, buffer);
  buffer=gsl_vector_ptr(state->dfx[i_prior],0);
  h5err&=H5LTmake_dataset_double(file_id,"LogPriorGradient", 1, size, buffer);
  // fisher information matrices
  size[0]=state->Hfx[i_posterior]->size1;
  size[1]=state->Hfx[i_posterior]->size2;
  buffer=gsl_matrix_ptr(state->Hfx[i_posterior],0,0);
  h5err&=H5LTmake_dataset_double(file_id,"PosteriorFisherInformation", 2, size, buffer);
  buffer=gsl_matrix_ptr(state->Hfx[i_likelihood],0,0);
  h5err&=H5LTmake_dataset_double(file_id,"LikelihoodFisherInformation", 2, size, buffer);
  buffer=gsl_matrix_ptr(state->Hfx[i_prior],0,0);
  h5err&=H5LTmake_dataset_double(file_id,"PriorFisherInformation", 2, size, buffer);
  // cholesky factor of the fisher information
  buffer=gsl_matrix_ptr(state->cholH_mat,0,0);
  h5err&=H5LTmake_dataset_double(file_id,"CholeskyFactorPostFI", 2, size, buffer);
  // close file
  h5err&=H5Fclose(file_id);
  assert(h5err==0);
  return EXIT_SUCCESS; 
}

int load_resume_state(const char *file_name, int rank, int R, const mcmc_kernel *kernel){
  assert(rank<R);
  assert(kernel);
  assert(file_name);
  hid_t file_id = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
  herr_t h5err=0;
  double *buffer;
  int original_rank=rank, original_comm_size=R;
  smmala_params *state = (smmala_params*) kernel->kernel_params;
  assert(state);
  // beta; the relaxation parameter of parallel tempering
  buffer=&(state->beta);
  h5err&=H5LTread_dataset_double(file_id,"beta",buffer);
  // step size;
  buffer=&(state->stepsize);
  h5err&=H5LTread_dataset_double(file_id,"stepsize",buffer);
  // MPI rank
  h5err&=H5LTread_dataset_int(file_id,"MPI_Comm_rank", &original_rank);
  h5err&=H5LTget_attribute_int(file_id,"MPI_Comm_rank", "MPI_Comm_size", &original_comm_size);
  assert(rank==original_rank);
  assert(R==original_comm_size);
  // posterior value fx
  // (likelihood, prior)
  buffer=state->fx;
  h5err&=H5LTread_dataset_double(file_id,"Posterior",buffer);
  assert(h5err==0);
  // state of the Markov chain x
  buffer=gsl_vector_ptr(state->x,0);
  h5err&=H5LTread_dataset_double(file_id,"MarkovChainState",buffer);
  // Gradients 
  buffer=gsl_vector_ptr(state->dfx[i_posterior],0);
  h5err&=H5LTread_dataset_double(file_id,"LogPosteriorGradient",buffer);
  buffer=gsl_vector_ptr(state->dfx[i_likelihood],0);
  h5err&=H5LTread_dataset_double(file_id,"LogLikelihoodGradient",buffer);
  buffer=gsl_vector_ptr(state->dfx[i_prior],0);
  h5err&=H5LTread_dataset_double(file_id,"LogPriorGradient",buffer);
  // fisher information matrices
  buffer=gsl_matrix_ptr(state->Hfx[i_posterior],0,0);
  h5err&=H5LTread_dataset_double(file_id,"PosteriorFisherInformation",buffer);
  buffer=gsl_matrix_ptr(state->Hfx[i_likelihood],0,0);
  h5err&=H5LTread_dataset_double(file_id,"LikelihoodFisherInformation",buffer);
  buffer=gsl_matrix_ptr(state->Hfx[i_prior],0,0);
  h5err&=H5LTread_dataset_double(file_id,"PriorFisherInformation",buffer);
  // cholesky factor of the fisher information
  buffer=gsl_matrix_ptr(state->cholH_mat,0,0);
  h5err&=H5LTread_dataset_double(file_id,"CholeskyFactorPostFI",buffer);
  // close file
  h5err&=H5Fclose(file_id);
  assert(h5err==0);
  return EXIT_SUCCESS; 
}



static smmala_params* smmala_params_alloc(const double beta, const int N, double step_size, const double target_acceptance){
  int overall_error=EXIT_SUCCESS;
  smmala_params* params = (smmala_params*) malloc(sizeof(smmala_params));
  if (params == NULL){
    fprintf(stderr,"smmala_params_alloc: malloc failed to allocate memory for smmala_params \n");
    return NULL;
  }
  /* There are three components to most values: an overall posterior value, a likelihood contribution to this value and a prior contribution
   * posterior(x) = beta*likelihood(x) + prior(x)
   * fx[0] = beta*fx[1] + fx[2]
   * nc=3
   */
  int i,nc=3;
  params->beta=beta;
  params->target_acceptance=target_acceptance;
  params->x=NULL;
  params->fx=NULL;
  params->dfx = NULL;
  params->Hfx = NULL;
  params->new_x = NULL;
  params->new_dfx = NULL;
  params->new_Hfx = NULL;
  params->mean_vec = NULL;
  params->cholH_mat = NULL;
  //printf("# [smmala_params_alloc] allocating memory.\n");
  params->fx=malloc(nc*sizeof(double));
  params->x = gsl_vector_alloc(N);
  if (params->x == NULL) overall_error++;
  
  params->dfx = malloc(nc*sizeof(gsl_vector*));
  for (i=0;i<nc;i++) params->dfx[i]=gsl_vector_alloc(N);
  if (params->dfx == NULL) overall_error++;

  params->Hfx = malloc(nc*sizeof(gsl_matrix*));
  for (i=0;i<nc;i++) params->Hfx[i]=gsl_matrix_alloc(N,N);
  if (params->Hfx == NULL) overall_error++;

  params->new_x = gsl_vector_alloc(N);
  if (params->new_x == NULL) overall_error++;
  
  params->new_dfx = malloc(nc*sizeof(gsl_vector*));
  for (i=0;i<nc;i++) params->new_dfx[i]=gsl_vector_alloc(N);
  if (params->new_dfx == NULL) overall_error++;
	
  params->new_Hfx = malloc(nc*sizeof(gsl_matrix*));
  for (i=0;i<nc;i++) params->new_Hfx[i]=gsl_matrix_alloc(N,N);
  if (params->new_Hfx == NULL) overall_error++;

  params->mean_vec = gsl_vector_alloc(N);
  if (params->mean_vec == NULL) overall_error++;
	
  params->cholH_mat = gsl_matrix_alloc(N,N);
  if (params->cholH_mat == NULL) overall_error++;

  params->fx[0] = 0;
  params->stepsize = step_size;
  //printf("allocations done: overall error: %i\n",overall_error);
  //fflush(stdout);
  //fflush(stderr);
  if (overall_error!=EXIT_SUCCESS) {
	fprintf(stderr,"[smmala_params_alloc] malloc failed \n");
	if (params->dfx) gsl_vector_free(params->x);	
	if (params->dfx){
	  for (i=0;i<nc;i++) if (params->dfx[i]!=NULL) gsl_vector_free(params->dfx[i]);
	  free(params->dfx);
	}
	if (params->Hfx){
	  for (i=0;i<nc;i++) if (params->Hfx[i]!=NULL) gsl_matrix_free(params->Hfx[i]);	  
	  free(params->Hfx);
	}
	
	if (params->new_x) gsl_vector_free(params->new_x);
	if (params->new_dfx) {
	  for (i=0;i<nc;i++) if (params->new_dfx[i]!=NULL) gsl_vector_free(params->new_dfx[i]);
	  free(params->new_dfx);
	}
	if (params->new_Hfx){
	  for (i=0;i<nc;i++) if(params->new_Hfx[i]!=NULL) gsl_matrix_free(params->new_Hfx[i]);
	  free(params->new_Hfx);
	}
	
	if (params->mean_vec) gsl_vector_free(params->mean_vec);
	if (params->cholH_mat) gsl_matrix_free(params->cholH_mat);
	free(params);

	return NULL;	
  }
  // printf("# [smmala_params_alloc] done.\n");
  //fflush(stdout);
  return(params);
}

static void smmala_params_free(smmala_params* params){
  int i,nc=3; // nc is the number of components to fx, dfx and Hfx;
  for (i=0;i<nc;i++) {
    gsl_vector_free(params->dfx[i]);
    gsl_vector_free(params->new_dfx[i]);
    gsl_matrix_free(params->Hfx[i]);
    gsl_matrix_free(params->new_Hfx[i]);
  }
  gsl_vector_free(params->x);
  gsl_vector_free(params->new_x);
  gsl_vector_free(params->mean_vec);
  gsl_matrix_free(params->cholH_mat);
  free(params->new_Hfx);
  free(params->new_dfx);
  free(params->Hfx);
  free(params->dfx);
  free(params->fx);
  free(params);
}

static int smmala_kernel_init(mcmc_kernel* kernel, const double *x){
  int res,i,n;
  n = kernel->N;
	
  smmala_params* state = (smmala_params*) kernel->kernel_params;
  /* copy x to the kernel x state */
  for (i=0; i<n; i++) kernel->x[i] = x[i];

  smmala_model* model = kernel->model_function;
  gsl_vector_const_view x_init = gsl_vector_const_view_array(x,n);
  res = model->LogPosterior(state->beta, &(x_init.vector), model->m_params, state->fx, state->dfx, state->Hfx);
  if (res != 0){
    fprintf(stderr,"smmala_kernel_init: Likelihood function failed\n");
    exit(-1);
  }
  res = gsl_linalg_cholesky_decomp(state->Hfx[0]);
  if (res != 0){
	fprintf(stderr,"Error: matrix not positive definite in smmala_init.\n");
	return -1;
  }  
  return EXIT_SUCCESS;
}

static int smmala_kernel_init_rand(mcmc_kernel* kernel){
  int res,i,n;
  n = kernel->N;
	
  gsl_rng* rng = (gsl_rng*) kernel->rng;
  smmala_model* model = kernel->model_function;
	
  /* sample random x from the prior */
  gsl_vector *x;
  x=gsl_vector_alloc(n);		/* automatic alloc */
  model->Prior_rnd(rng, model->m_params, x->data);
  
  /* copy x to the kernel x state */
  for ( i=0; i < n; i++)
    kernel->x[i] = gsl_vector_get(x,i);
  
  smmala_params* state = (smmala_params*) kernel->kernel_params;
  
  res = model->LogPosterior(state->beta, x, model->m_params, state->fx,state->dfx, state->Hfx);
  /* TODO: write a proper error handler */
  if (res != 0){
    fprintf(stderr,"smmala_kernel_init: Likelihood function failed\n");
    return 1;
  }
  
  res = gsl_linalg_cholesky_decomp(state->Hfx[0]);
  if (res != 0){
    fprintf(stderr,"Error: matrix not positive definite in smmala_init.\n");
    return -1;
  }
  return 0;
}


static int nat_grad_step(const gsl_vector* x, const gsl_matrix* cholPr,  const gsl_vector* grad, gsl_vector* result, double stepSize){
  /* newx = x + (0.5*e)*H^{-1}*grad */
  /* H^{-1}*grad is similar to solving the system H*result = grad */
  gsl_linalg_cholesky_solve(cholPr, grad, result);
  gsl_vector_scale(result, 0.5*stepSize);
  gsl_vector_add(result, x);
  
  return 0;
}


static int smmala_kernel_sample(mcmc_kernel* kernel, int* acc){
  int i;
  smmala_params* state = (smmala_params*) kernel->kernel_params;
  smmala_model* model = (smmala_model*) kernel->model_function;
  gsl_rng* rng = (gsl_rng*) kernel->rng;
  
  //int n = kernel->N;
  double stepsize = state->stepsize;
  
  gsl_matrix_memcpy(state->cholH_mat, state->Hfx[0]);
  
  /* Calculate mean vector for random (Gaussian) step, by following the current gradient of posterior */
  nat_grad_step(state->x, state->Hfx[0], state->dfx[0], state->mean_vec,stepsize);
  /* Random vector from multivariate normal with mean mean_vec and precision cholH_mat */
  gsl_matrix_scale(state->cholH_mat,1.0/sqrt(stepsize));
  mnv_norm_rnd_cholPr(rng, state->mean_vec, state->cholH_mat, state->new_x);
  /* Calculate propbability of new state given old state */
  double pNgO = log_mv_norm_pdf_cholP(state->new_x, state->mean_vec, state->cholH_mat);
  
  /* evaluate model and get new state */
  double new_fx[3];
  int res = model->LogPosterior(state->beta,
				state->new_x,
				model->m_params,
				new_fx,
				state->new_dfx,
				state->new_Hfx);
  if (res != GSL_SUCCESS){
    fprintf(stderr,"[smmala warning]: Model cannot be evaluated with this argument: Likelihood=0\n");
    *acc = 0;
  } else {
    res = gsl_linalg_cholesky_decomp(state->new_Hfx[0]);
    if (res != 0){
	fprintf(stderr,"Warning: matrix not positive definite in smmala_sample.\n Will replace matrix with I (gsl_matrix_set_identity).\n");
	gsl_matrix_set_identity(state->new_Hfx[0]);
    }
    gsl_matrix_memcpy(state->cholH_mat, state->new_Hfx[0]);	
    /* Calculate new mean */
    nat_grad_step(state->new_x, state->cholH_mat, state->new_dfx[0], state->mean_vec, stepsize);
    
    /* Calculate propbability of old state given new state */
    gsl_matrix_scale(state->cholH_mat, 1.0/sqrt(stepsize));
    double pOgN = log_mv_norm_pdf_cholP(state->x, state->mean_vec, state->cholH_mat);
    
    /* Accept/Reject new state */
    double mh_ratio = new_fx[0] + pOgN - state->fx[0] - pNgO;
    double rand_dec = log(gsl_rng_uniform(rng));
    if ( (mh_ratio > 0.0)||(mh_ratio > rand_dec) ) {
      *acc = 1;
      gsl_vector_memcpy(state->x,state->new_x);
      memcpy(state->fx, new_fx,sizeof(double)*3);
      for (i=0;i<3;i++){
	gsl_vector_memcpy(state->dfx[i],state->new_dfx[i]);
	gsl_matrix_memcpy(state->Hfx[i],state->new_Hfx[i]);
      }
    }else{
      *acc = 0;
    }
  }
  return GSL_SUCCESS;
}

static void smmala_kernel_adapt(mcmc_kernel* kernel, double acc_rate){
  smmala_params* params = (smmala_params*) kernel->kernel_params;
  double a=params->target_acceptance;
  if (acc_rate > 1.1*a) {
    params->stepsize = 1.3*params->stepsize;
  }
  else if(acc_rate < 0.9*a){
    params->stepsize = 0.8*params->stepsize;
  }  
}

static void smmala_kernel_free(mcmc_kernel* kernel){
  smmala_params* params = (smmala_params*)kernel->kernel_params;
  smmala_params_free(params);
  gsl_rng* rng = (gsl_rng*)kernel->rng;
  gsl_rng_free(rng);
  free(kernel);
}

static void smmala_kernel_print_stats(mcmc_kernel* kernel, FILE* outStream){
  smmala_params* params = (smmala_params*)kernel->kernel_params;
  fprintf(outStream, "Step-size: %5.4g\n",params->stepsize);
}

smmala_model* smmala_model_alloc(fptrPosterior_smmala Lx, fptrPrior_rnd Prx,  void* model_params){
  smmala_model* model = (smmala_model*) malloc(sizeof(smmala_model));
  if (model == NULL){
    fprintf(stderr,"malloc failed to allocate memory for smmala_model\n");
    exit(-1);
  }
  model->LogPosterior = Lx;
  model->Prior_rnd = Prx;
  model->m_params = model_params;
  return model;
}

void smmala_model_free(smmala_model* model){
  free(model);
}

double smmala_kernel_get_beta(mcmc_kernel *kernel){
  smmala_params *state=kernel->kernel_params;
  return state->beta;
}

mcmc_kernel* smmala_kernel_alloc(const double beta, const int N, double step_size, smmala_model* model_function, const unsigned long int seed, const double target_acceptance){
  smmala_params* state = smmala_params_alloc(beta, N, step_size, target_acceptance);
  if( state == NULL ){
    fprintf(stderr,"malloc failed to allocate memory for params in smmala_kernel_alloc \n");
    exit(-1);
  }
  
  mcmc_kernel* kernel = (mcmc_kernel*) malloc(sizeof(mcmc_kernel));
  if (kernel == NULL){
    smmala_params_free(state);
    fprintf(stderr,"malloc failed to allocate memory for mcmc_kernel in smmala_alloc.\n");
    exit(-1);
  }

  kernel->x=NULL;
  kernel->fx=NULL;

  if (state->x->data!=NULL) kernel->x = state->x->data;
  if (state->fx!=NULL) kernel->fx = state->fx;
  
  if (kernel->x == NULL || kernel->fx == NULL){
    smmala_params_free(state);
    free(kernel);
    fprintf(stderr,"smmala: kernel->x OR kernel->fx is NULL pointer, check smmala parameter allocation \n");
    exit(-1);
  }
  
  gsl_rng_env_setup();
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  if (r == NULL){
    smmala_params_free(state);
    free(kernel);
    fprintf(stderr,"failed to create gsl_rng in smmala_alloc \n");
    exit(-1);
  }
  
  gsl_rng_set(r, seed);
  kernel->rng = r;
  
  kernel->N = N;
  kernel->model_function = model_function;
  kernel->kernel_params = state;
  kernel->ExchangeInformation = &smmala_exchange_information;
  kernel->SwapChains = &smmala_swap_chains;
  kernel->Sample = &smmala_kernel_sample;
  kernel->Adapt = &smmala_kernel_adapt;
  kernel->Init = &smmala_kernel_init;
  kernel->InitR = &smmala_kernel_init_rand;
  kernel->Free = &smmala_kernel_free;
  kernel->PrintStats = &smmala_kernel_print_stats;
  kernel->GetBeta = &smmala_kernel_get_beta;
  return kernel;
}
