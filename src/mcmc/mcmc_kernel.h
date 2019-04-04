/*
 *  mcmc_kernel.h
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

/*
 * MCMC kernel interface
 * Complie with -DHAVE_INLINE flag to enable inline optimisations used in this module.
 */

#ifndef __MCMC_KERNEL_H__
#define __MCMC_KERNEL_H__

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>

  /* Abstract data type for all MCMC kernels. All MCMC algorithms implement this interface */
  typedef struct mcmc_kernel_struct_ mcmc_kernel;

  /* ====== Common interface for all MCMC samplers ======== */

  /* Draw one sample from the Markov Chain.
   * Arguments:
   *	kernel:	    a pointer to the MCMC kernel data structure.
   *  Result:
   *	 returns zero for success and non-zero for failure.
   *	 the new sample is directly updated in kernel->x.
   *	 acc is set to 0 if the chain  in the previous state (reject)
   *   and 1 if the chain made a transition to a new state (accept)
   */
  typedef	int (*fptrMCMCSample)(mcmc_kernel* kernel, int* acc);
  int mcmc_sample(mcmc_kernel* kernel, int* acc);

  /* This function is used for Parallel Tempering communication.
   *
   */
  typedef int (*fptrMCMCExchangeInformation)(mcmc_kernel* kernel, const int DEST, void *buffers);
  int mcmc_exchange_information(mcmc_kernel* kernel, const int DEST, void *buffer);
  //smmala_exchange_information(mcmc_kernel* kernel, const int DEST, void *buffer)

  typedef int (*fptrMCMCSwapChains)(mcmc_kernel* kernel, const int master, const int rank, const int DEST, void *buffer);
  int mcmc_swap_chains(mcmc_kernel* kernel, const int master, const int rank, const int DEST, void *buffer);
  //smmala_swap_chains(mcmc_kernel* kernel, const int master, const int rank, const int other_rank, void *buffer)
  
  /* Initialise MCMC kernel with initial parameters x.
   * Arguments:
   *	kernel:		a pointer to the MCMC kernel structure.
   *	x:			an array of N doubles with initial parameters. N must be 
   *				equal to kernel->N.
   * Result:
   * returns 0 for success and non-zero for error.
   */	
  typedef	int (*fptrMCMCInit)(mcmc_kernel* kernel, const double* x);
  int mcmc_init(mcmc_kernel* kernel, const double* x);

  /* Initialise MCMC kernel using random initial parameters sampled from the prior.
   * Arguments:
   *	kernel:		a pointer to the MCMC kernel structure.
   * Result:
   *	returns 0 for success and non-zero for error.
   */	
  typedef	int (*fptrMCMCInitRand)(mcmc_kernel* kernel);
  int mcmc_init_rand(mcmc_kernel* kernel);	
	
  /* Free MCMC data structure and de-allocate all memory.
   * Arguments:
   *	kernel:		a pointer to the MCMC kernel structure.
   */
  typedef	void (*fptrMCMCFree)(mcmc_kernel* kernel);
  void mcmc_free(mcmc_kernel* kernel);
	
  /* Adapt MCMC parameters during the burn-in.
   * Arguments:
   *	kernel:		a pointer to the MCMC kernel structure.
   *  acc_rate:   acceptance rate.
   */
  typedef	void (*fptrMCMCAdapt)(mcmc_kernel* kernel, double acc_rate);
  void mcmc_adapt(mcmc_kernel* kernel, double acc_rate);
	
  /* Print statistics and information for the MCMC sampler.
   * Arguments:
   *	kernel:		a pointer to the MCMC kernel structure.
   *  outStream:  an output stream to write.
   */	
  typedef	void (*fptrMCMCPrint)(mcmc_kernel* kernel, FILE* outStream);
  void mcmc_print_stats(mcmc_kernel* kernel, FILE* s);

  /* Common function for all mcmc kernles that prints the current sample to an output stream 
   * Arguments:
   *	kernel:		a pointer to the MCMC kernel structure.
   *  outStream:  an output stream to write.
   */
  void mcmc_print_sample(mcmc_kernel* kernel, FILE* s);
  int mcmc_write_sample(mcmc_kernel *kernel, FILE *s);
  typedef double (*fptrGetBeta)(mcmc_kernel* kernel);
  double mcmc_get_beta(mcmc_kernel *kernel);
  /* ==================================================================== */

  /* Implementation details */
  /* TODO: hide implementation details from public header file */

  /* Data structure for every mcmc kernel */
  struct mcmc_kernel_struct_ {
    int N;						/* number of parameters. size of x */
    long int t;					/* itteration number, initialise at 0 */
    double* x;					/* link to method specific state, e.g.: sth. like smmala->state->x->data; current state of the Markov chain, ~ ODE-model log-parameters */
    double* fx;                                 /* this is also a link to the LogPosterior value saved in the working memory of the specific mcmc method */
    void* rng;					/* pointer to random number generator */
    void* model_function;		/* pointer to structure with user defined functions for the model */
    void* kernel_params;		/* kernel specific parameters and working storage */
    fptrMCMCExchangeInformation ExchangeInformation;
    fptrMCMCSwapChains SwapChains;
    fptrMCMCSample Sample;		/* sampling function for mcmc kernel */
    fptrMCMCInit Init;			/* initialisation function for mcmc kernel using predifined starting values for parameters */
    fptrMCMCInitRand InitR;		/* initialisation function for mcmc kernel using random sample from the prior */
    fptrMCMCFree Free;			/* de-allocation function for mcmc kernel. Frees all memory used by the mcmc kernel */ 
    fptrMCMCAdapt Adapt;		/* adaptation function for mcmc kernel. Addapts any tuning kernel parameters */
    fptrMCMCPrint PrintStats;	/* loging function for mcmc kernel. Prints to an output stream kernel statistics and loging information */
    fptrGetBeta GetBeta;
  };
  
  /* Draw a radom sample from the prior.
   *  Arguments:
   *   rng:			a pointer to a GSL random number generator.
   *	 model_params:  a pointer to a user defined model parameters structure.
   *  Results:
   *	 sample is the array with a random sample from the prior.
   */
  typedef int (*fptrPrior_rnd) (gsl_rng* rng, const void* model_params, double* sample);
  
  /* convinience macros for sampling and adaptation */
#define MCMC_SAMPLE( K, a ) ((K)->Sample) ((K), (a))
#define MCMC_ADAPT( K, a ) ((K)->Adapt) ((K), (a))
#define MCMC_DIM(K) (K->N)
#define MCMC_STATE(K) (K->x)
#define MCMC_POSTERIOR(K) (K->fx)
  
#ifdef __cplusplus
}
#endif

#endif
