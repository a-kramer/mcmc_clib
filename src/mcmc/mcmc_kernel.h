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
 */

#ifndef _MCMC_KERNEL_H
#define _MCMC_KERNEL_H

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>

#define MCMC_SUCCESS 0
#define MCMC_POSTERIOR_FAILED -1
#define MCMC_METHOD_FAILED -2

typedef struct mcmc_kernel_struct_ mcmc_kernel;
int mcmc_sample(mcmc_kernel* kernel, int* acc);
int mcmc_exchange_information(mcmc_kernel* kernel, const int DEST, void *buffer);
buffer);
int mcmc_swap_chains(mcmc_kernel* kernel, const int master, const int rank, const int DEST, void *buffer);
int mcmc_init(mcmc_kernel* kernel, const double* x);
int mcmc_init_rand(mcmc_kernel* kernel);	
mcmc_kernel* mcmc_alloc(const double beta, const int N, double step_size, statistical_model* smod, const unsigned long int seed, const double target_acceptance);
void mcmc_free(mcmc_kernel* kernel);
void mcmc_adapt(mcmc_kernel* kernel, double acc_rate);
void mcmc_print_stats(mcmc_kernel* kernel, FILE* s);
void mcmc_print_sample(mcmc_kernel* kernel, FILE* s);
int mcmc_write_sample(mcmc_kernel *kernel, FILE *s);
struct mcmc_kernel_struct_ {
  long int iter; 
  double beta;
  gsl_vector* x;
  double* fx;
  void* rng;
  void* statistical_model; 
  void* mcmc_specific_params;
};
#endif
