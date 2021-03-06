/*
 *  mcmc_kernel.c
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
#include <stdlib.h>
#include "mcmc_kernel.h"
#include "hdf5.h"
#include "hdf5_hl.h"


void mcmc_print_sample(mcmc_kernel* kernel, FILE* s){
  int n = kernel->N;
  int i;
  if (kernel->x!=NULL && kernel->fx!=NULL){
    for (i = 0; i < n; i++) {
      fprintf(s, "%12.5e ", kernel->x[i]);
    }
    fprintf(s," %12.5e ",kernel->fx[0]);
    fprintf(s, "\n");
  } else {
    perror("kernel->(x or fx) is a NULL pointer.");
    exit(-1);
  }
}

int mcmc_write_sample(mcmc_kernel *kernel, FILE *s){
  int n = kernel->N;
  fwrite(kernel->x,sizeof(double),n,s);
  fwrite(kernel->fx,sizeof(double),1,s);
  return EXIT_SUCCESS;
}


int mcmc_exchange_information(mcmc_kernel* kernel, const int DEST, void *buffer){
  return kernel->ExchangeInformation(kernel, DEST, buffer);
}

int mcmc_swap_chains(mcmc_kernel* kernel, const int master, const int rank, const int DEST, void *buffer){
  return kernel->SwapChains(kernel, master, rank, DEST, buffer);
}
  

/* used only when the code is not compiled with the -DHAVE_INLINE flag */
int mcmc_sample(mcmc_kernel* kernel, int* acc){
	kernel->t ++;
	return kernel->Sample(kernel, acc);
}

void mcmc_adapt(mcmc_kernel* kernel, double acc_rate){
	kernel->Adapt(kernel, acc_rate);
}

int mcmc_init(mcmc_kernel* kernel, const double *x){
	return kernel->Init(kernel, x);
}

int mcmc_init_rand(mcmc_kernel* kernel){
	return kernel->InitR(kernel);
}

void mcmc_free(mcmc_kernel* kernel){
	kernel->Free(kernel);
}

void mcmc_print_stats(mcmc_kernel* kernel, FILE* s){
	kernel->PrintStats(kernel, s);
}

double mcmc_get_beta(mcmc_kernel *kernel){
  return kernel->GetBeta(kernel);
}
