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
#include "mcmc_kernel.h"
#include "mv_norm.h"
#include "model_parameters_smmala.h"
#include "../app/diagnosis_output.h"
#include "hdf5.h"
#include "hdf5_hl.h"


#define SWAP(a, b, tmp)  tmp = a; a = b; b = tmp

// fx = beta*lx + px
typedef struct {
  double stepsize;
  gsl_vector **dfx;	/* log-posterior gradient wrt to x*/
  gsl_matrix **Hfx;     /* cholesky factor of Fisher Information*/
  gsl_vector *new_x;	/* tmp working space for new x */
  gsl_vector **new_dfx;
  gsl_matrix **new_Hfx;
  gsl_vector *mean_vec;
  gsl_matrix *cholH_mat;
  double target_acceptance;
} smmala_params;

static smmala_params* smmala_params_alloc(const double beta, const int N, double step_size, const double target_acceptance);

double get_step_size(const mcmc_kernel *kernel){
  smmala_params *state = kernel->mcmc_specific_params;
  return state->stepsize;
}

static void print_on_error(int error_code, const char *attempted_comm){
  char error_string[MPI_MAX_ERROR_STRING];
  int len;
  if (error_code != MPI_SUCCESS){
    MPI_Error_string(ec, error_string, &len);
    fprintf(stderr, "%s, error: %s\n", attempted_comm, error_string);
    MPI_Error_class(ec, &error_class);
    MPI_Error_string(error_class, error_string, &len);
    fprintf(stderr, "%s, class: %s\n", attempted_comm, error_string);
  }
}

static void mcmc_exchange_information(mcmc_kernel* kernel, const int DEST, double *buffer, size_t len){
  MPI_Status status; // MPI_Status contains: MPI_SOURCE, MPI_TAG, MPI_ERROR
  int N=kernel->x->size; // size of x, dfx; N*N is size of FI;
  smmala_params *state = (smmala_params*) kernel->mcmc_specific_params;
  int TAG=0;
  int SRC=DEST; /* send to and receive from the same process */

  int size=1;
  int ec;
  assert(len>=2);
  ec=MPI_Sendrecv(&(kernel->beta), size, MPI_DOUBLE, DEST, TAG,
		  buffer[0], size, MPI_DOUBLE, SRC, TAG,
		  MPI_COMM_WORLD, &status);
  print_on_error(ec,"sendrecv beta");

  ec=MPI_Sendrecv(&(kernel->fx[1]), size, MPI_DOUBLE, DEST, TAG,
		  buffer[1], size, MPI_DOUBLE, SRC, TAG,
		  MPI_COMM_WORLD, &status);
  print_on_error(ec,"sendrecv likelihood");
}

int mcmc_swap_chains(mcmc_kernel* kernel, const int is_sender, const int rank, const int other_rank){
  double a;
  int swap_accepted=0;
  //int N=kernel->N; // size of x, dfx; N*N is size of FI;
  gsl_rng *rng = (gsl_rng*) kernel->rng;
  smmala_params *state = (smmala_params*) kernel->mcmc_specific_params;
  double buffer[2];
  mcmc_exchange_information(kernel, other_rank, buffer, 2);
  int TAG=0;
  double beta=kernel->beta;
  
  assert(len>=2);
  double their_beta=buffer[0];
  double their_fx_1=buffer[1];

  int gsl_status;
  gsl_sf_result result;
  double beta_state = (their_beta - beta)*(state->fx[1] - their_fx_1);
  
  gsl_status=gsl_sf_exp_e(beta_state, &result);
  if (gsl_status==GSL_SUCCESS && result.err<result.val){
    a=result.val;
  } else {
    a=0.0;
  }
  
  if (is_sender){
    double r1=gsl_rng_uniform(rng);
    swap_accepted=(r1<a);
    MPI_Send(&swap_accepted, 1, MPI_INT, other_rank,TAG,MPI_COMM_WORLD);
  }else{
    MPI_Recv(&swap_accepted, 1, MPI_INT, other_rank,TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (swap_accepted){
    /* swap the two beta values */
    beta=their_beta;
    kernel->beta=beta;
    /* re-calculate beta-dependent values */
    state->fx[0]=beta*state->fx[1] + state->fx[2];
    gsl_vector_memcpy(state->dfx[0],state->dfx[1]);
    gsl_vector_scale(state->dfx[0],beta);
    gsl_vector_add(state->dfx[0],state->dfx[2]);

    gsl_matrix_memcpy(state->Hfx[0],state->Hfx[1]);
    gsl_matrix_scale(state->Hfx[0],beta*beta);
    gsl_matrix_add(state->Hfx[0],state->Hfx[2]);
    assert(gsl_linalg_cholesky_decomp(state->Hfx[0])==GSL_SUCCESS);
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
  hsize_t size[2];
  herr_t h5err=0;
  double *buffer;
  smmala_params *state = kernel->mcmc_specific_params;
  assert(state);
  
  size[0]=1;
  buffer=&(state->beta);
  h5err|=H5LTmake_dataset_double(file_id,"beta", 1, size, buffer);

  buffer=&(state->stepsize);
  h5err|=H5LTmake_dataset_double(file_id,"stepsize", 1, size, buffer);
  
  h5err|=H5LTmake_dataset_int(file_id,"MPI_Comm_rank", 1, size, &rank);
  h5err|=H5LTset_attribute_int(file_id,"MPI_Comm_rank", "MPI_Comm_size", &R,1);
  
  h5err|=H5LTset_dataset_string(file_id,"mcmc method", "SMMALA");
  
  size[0]=3;
  buffer=kernel->fx;
  h5err|=H5LTmake_dataset_double(file_id,"Posterior", 1, size, buffer);
  h5err|=H5LTset_attribute_string(file_id,"Posterior", "key", "LogPosterior; LogLikelihood; LogPrior");
  h5err|=H5LTset_attribute_string(file_id,"Posterior", "info", "LogPosterior=beta*LogLikelihood+LogPrior");
  assert(h5err==0);

  // state of the Markov chain x
  size[0]=kernel->x->size;
  buffer=gsl_vector_ptr(kernel->x,0);
  h5err|=H5LTmake_dataset_double(file_id,"MarkovChainState", 1, size, buffer);
  /* gradients */
  size[0]=state->dfx[i_posterior]->size;
  buffer=gsl_vector_ptr(state->dfx[i_posterior],0);
  h5err|=H5LTmake_dataset_double(file_id,"LogPosteriorGradient", 1, size, buffer);
  buffer=gsl_vector_ptr(state->dfx[i_likelihood],0);
  h5err|=H5LTmake_dataset_double(file_id,"LogLikelihoodGradient", 1, size, buffer);
  buffer=gsl_vector_ptr(state->dfx[i_prior],0);
  h5err|=H5LTmake_dataset_double(file_id,"LogPriorGradient", 1, size, buffer);
  /* fisher information matrices */
  size[0]=state->Hfx[i_posterior]->size1;
  size[1]=state->Hfx[i_posterior]->size2;
  buffer=gsl_matrix_ptr(state->Hfx[i_posterior],0,0);
  h5err|=H5LTmake_dataset_double(file_id,"PosteriorFisherInformation", 2, size, buffer);
  buffer=gsl_matrix_ptr(state->Hfx[i_likelihood],0,0);
  h5err|=H5LTmake_dataset_double(file_id,"LikelihoodFisherInformation", 2, size, buffer);
  buffer=gsl_matrix_ptr(state->Hfx[i_prior],0,0);
  h5err|=H5LTmake_dataset_double(file_id,"PriorFisherInformation", 2, size, buffer);
  /* cholesky factor of the fisher information */
  buffer=gsl_matrix_ptr(state->cholH_mat,0,0);
  h5err|=H5LTmake_dataset_double(file_id,"CholeskyFactorPostFI", 2, size, buffer);

  h5err|=H5Fclose(file_id);
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
  buffer=&(kernel->beta);
  h5err|=H5LTread_dataset_double(file_id,"beta",buffer);
  // step size;
  buffer=&(state->stepsize);
  h5err|=H5LTread_dataset_double(file_id,"stepsize",buffer);
  // MPI rank
  h5err|=H5LTread_dataset_int(file_id,"MPI_Comm_rank", &original_rank);
  h5err|=H5LTget_attribute_int(file_id,"MPI_Comm_rank", "MPI_Comm_size", &original_comm_size);
  assert(rank==original_rank);
  assert(R==original_comm_size);
  // posterior value fx
  // (likelihood, prior)
  buffer=state->fx;
  h5err|=H5LTread_dataset_double(file_id,"Posterior",buffer);
  assert(h5err==0);
  // state of the Markov chain x
  buffer=gsl_vector_ptr(state->x,0);
  h5err|=H5LTread_dataset_double(file_id,"MarkovChainState",buffer);
  // Gradients 
  buffer=gsl_vector_ptr(state->dfx[i_posterior],0);
  h5err|=H5LTread_dataset_double(file_id,"LogPosteriorGradient",buffer);
  buffer=gsl_vector_ptr(state->dfx[i_likelihood],0);
  h5err|=H5LTread_dataset_double(file_id,"LogLikelihoodGradient",buffer);
  buffer=gsl_vector_ptr(state->dfx[i_prior],0);
  h5err|=H5LTread_dataset_double(file_id,"LogPriorGradient",buffer);
  // fisher information matrices
  buffer=gsl_matrix_ptr(state->Hfx[i_posterior],0,0);
  h5err|=H5LTread_dataset_double(file_id,"PosteriorFisherInformation",buffer);
  buffer=gsl_matrix_ptr(state->Hfx[i_likelihood],0,0);
  h5err|=H5LTread_dataset_double(file_id,"LikelihoodFisherInformation",buffer);
  buffer=gsl_matrix_ptr(state->Hfx[i_prior],0,0);
  h5err|=H5LTread_dataset_double(file_id,"PriorFisherInformation",buffer);
  // cholesky factor of the fisher information
  buffer=gsl_matrix_ptr(state->cholH_mat,0,0);
  h5err|=H5LTread_dataset_double(file_id,"CholeskyFactorPostFI",buffer);
  // close file
  h5err|=H5Fclose(file_id);
  assert(h5err==0);
  return EXIT_SUCCESS; 
}

/* There are three components to most values: an overall posterior
 * value, a likelihood contribution to this value and a prior
 * contribution posterior(x) = beta*likelihood(x) + prior(x) fx[0] =
 * beta*fx[1] + fx[2];
 */
static smmala_params* smmala_params_alloc(const double beta, const int N, double step_size, const double target_acceptance){
  int overall_error=EXIT_SUCCESS;
  smmala_params* params = (smmala_params*) malloc(sizeof(smmala_params));
  assert(params);
  int i,nc=3;
  params->target_acceptance=target_acceptance;
  params->dfx = malloc(nc*sizeof(gsl_vector*));
  for (i=0;i<nc;i++) params->dfx[i]=gsl_vector_alloc(N);

  params->Hfx = malloc(nc*sizeof(gsl_matrix*));
  for (i=0;i<nc;i++) params->Hfx[i]=gsl_matrix_alloc(N,N);

  params->new_x = gsl_vector_alloc(N);
  params->new_dfx = malloc(nc*sizeof(gsl_vector*));
  for (i=0;i<nc;i++) params->new_dfx[i]=gsl_vector_alloc(N);

  params->new_Hfx = malloc(nc*sizeof(gsl_matrix*));
  for (i=0;i<nc;i++) params->new_Hfx[i]=gsl_matrix_alloc(N,N);

  params->mean_vec = gsl_vector_alloc(N);
  params->cholH_mat = gsl_matrix_alloc(N,N);
  params->stepsize = step_size;
  return params;
}

void smmala_params_free(smmala_params* params){
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

int mcmc_kernel_init(mcmc_kernel* kernel, const gsl_vector *x){
  int res,i;
  assert(kernel);
  assert(x);
  smmala_params* state = (smmala_params*) kernel->mcmc_specific_params;
  gsl_vector_memcpy(kernel->x,x);
  statistical_model* model = kernel->statistical_model;
  assert(model);
  res = model->LogPosterior(state->beta, x, model->m_params, state->fx, state->dfx, state->Hfx);
  if (res){
    fprintf(stderr,"[%s] Likelihood function failed\n",__func__);
    abort();
  }
  res = gsl_linalg_cholesky_decomp(state->Hfx[0]);
  if (res){
    fprintf(stderr,"[%s] matrix not positive definite.\n",__func__);
    abort();
  }
  printf("[%s] done.\n",__func__);
  fflush(stdout);
  return EXIT_SUCCESS;
}

static int nat_grad_step(const gsl_vector* x, const gsl_matrix* cholPr,  const gsl_vector* grad, gsl_vector* result, double stepSize){
  /* newx = x + (0.5*e)*H^{-1}*grad */
  /* H^{-1}*grad is similar to solving the system H*result = grad */
  gsl_linalg_cholesky_solve(cholPr, grad, result);
  gsl_vector_scale(result, 0.5*stepSize);
  gsl_vector_add(result, x);
  return 0;
}


int mcmc_sample(mcmc_kernel* kernel, int* acc){
  int i;
  smmala_params *state = (smmala_params*) kernel->mcmc_specific_params;
  statistical_model *model = (statistical_model*) kernel->statistical_model;
  gsl_rng *rng = (gsl_rng*) kernel->rng;
  
  //int n = kernel->N;
  double stepsize = state->stepsize;
  
  gsl_matrix_memcpy(state->cholH_mat, state->Hfx[0]);
  
  /* Calculate mean vector for random (Gaussian) step, by following the current gradient of posterior */
  nat_grad_step(kernel->x, state->Hfx[0], state->dfx[0], state->mean_vec, stepsize);
  /* Random vector from multivariate normal with mean mean_vec and precision cholH_mat */
  gsl_matrix_scale(state->cholH_mat,1.0/sqrt(stepsize));
  mnv_norm_rnd_cholPr(rng, state->mean_vec, state->cholH_mat, state->new_x);
  /* Calculate propbability of new state given old state */
  double pNgO = log_mv_norm_pdf_cholP(state->new_x, state->mean_vec, state->cholH_mat);
  
  /* evaluate model and get new state */
  double new_fx[3];
  int res = model->LogPosterior
    (kernel->beta,
     state->new_x,
     model->m_params,
     new_fx,
     state->new_dfx,
     state->new_Hfx);

  if (res != GSL_SUCCESS){
    fprintf(stderr,"[smmala warning]: Model cannot be evaluated with this argument: Likelihood=0\n");
    *acc = 0;
    return MCMC_POSTERIOR_FAILURE;
  }
  
  res = gsl_linalg_cholesky_decomp(state->new_Hfx[0]);
  if (res){
    fprintf(stderr,"Warning: matrix not positive definite in smmala_sample.\n\
Calling gsl_matrix_set_identity().\n");
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
  if ((mh_ratio > 0.0) || (mh_ratio > rand_dec)) {
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
  return MCMC_SUCCESS;
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

static void mcmc_kernel_free(mcmc_kernel* kernel){
  smmala_params* params = (smmala_params*)kernel->kernel_params;
  smmala_params_free(params);
  gsl_rng* rng = (gsl_rng*)kernel->rng;
  gsl_rng_free(rng);
  free(kernel);
}

static void mcmc_print_stats(mcmc_kernel* kernel, FILE* outStream){
  smmala_params* p = (smmala_params*) kernel->mcmc_specific_params;
  fprintf(outStream, "Step-size: %5.4g\n",p->stepsize);
}

statistical_model* statistical_model_alloc(fptrPosterior_smmala Lx, fptrPrior_rnd Prx,  void* model_params){
  statistical_model* model = malloc(sizeof(smmala_model));
  assert(model);
  model->LogPosterior = Lx;
  model->Prior_rnd = Prx;
  model->m_params = model_params;
  return model;
}

mcmc_kernel* mcmc_kernel_alloc(const double beta, const int N, double step_size, statistical_model* smod, const unsigned long int seed, const double target_acceptance){
  smmala_params* state = smmala_params_alloc(beta, N, step_size, target_acceptance);
  mcmc_kernel* kernel = (mcmc_kernel*) malloc(sizeof(mcmc_kernel));
  assert(state);
  assert(kernel);
  kernel->iter=0;
  kernel->x=gsl_vector_alloc(N);
  kernel->fx=malloc(sizeof(double)*3);
  gsl_rng_env_setup();
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  assert(r);
  gsl_rng_set(r, seed);
  kernel->rng = r;
  kernel->statistical_model = smod;
  kernel->mcmc_specific_param = state;
  return kernel;
}
