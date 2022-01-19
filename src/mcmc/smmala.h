/*
 * SMMALA interface
 */

#ifndef __SMMALA_H__
#define __SMMALA_H__
  
#include <gsl/gsl_rng.h>
#include "mcmc_kernel.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define i_posterior 0
#define i_likelihood 1
#define i_prior 2
  
typedef int (*fptrPosterior_smmala)
(const double beta,
 const gsl_vector *x,
 void* model_params,
 double *fx,
 gsl_vector **dfx, gsl_matrix **Hfx);

/*typedef int (*fptrModelGrad)(const double* x, const void* model_params,
  double* dfx, double* Hfx);*/

typedef struct{
  fptrPosterior_smmala LogPosterior;
  fptrPrior_rnd  Prior_rnd;
  void* m_params;
} smmala_model;

void* smmala_comm_buffer_alloc(int D);  
smmala_model* smmala_model_alloc(fptrPosterior_smmala Lx, fptrPrior_rnd Prx, void* model_params);

void smmala_model_free(smmala_model* model);

mcmc_kernel* smmala_kernel_alloc(const double beta, const int N, double step_size, smmala_model* model_function, unsigned long int seed, const double target_acceptance);

int write_resume_state(const char *file_name, int rank, int R, const mcmc_kernel *kernel);
int load_resume_state(const char *file_name, int rank, int R, const mcmc_kernel *kernel);
double get_step_size(const mcmc_kernel *kernel);
#endif
