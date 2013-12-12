/*
 * SMMALA interface
 */

#ifndef __SMMALA_H__
#define __SMMALA_H__

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
	
#include <gsl/gsl_rng.h>
#include "mcmc_kernel.h"

  typedef int (*fptrPosterior_smmala)(const double* x, const void* model_params, double* fx,
								double* dfx, double* FIx);
	
  /*typedef int (*fptrModelGrad)(const double* x, const void* model_params,
							   double* dfx, double* Hfx);*/
	
  typedef struct{
	  fptrPosterior_smmala Likelihood;
	  fptrPrior_rnd  Prior_rnd;
	  void* m_params;
  } smmala_model;
		

  smmala_model* smmala_model_alloc(fptrPosterior_smmala Lx, fptrPrior_rnd Prx, void* model_params);

  void smmala_model_free(smmala_model* model);

  mcmc_kernel* smmala_kernel_alloc(int N, double step_size, smmala_model* model_function, unsigned long int seed, double target_acceptance);
  
  
#ifdef __cplusplus
}
#endif

#endif
