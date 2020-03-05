#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "../ode/ode_model.h"
#include "smmala.h"

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
  buffer=state->beta;
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
  buffer=gsl_matrix_ptr(state->cholH_mat[i_posterior],0,0);
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
  buffer=state->beta;
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
