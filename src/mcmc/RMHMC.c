/*
 *  RMHMC.c
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


#include <stdio.h>
#include "RMHMC.h"

/* Auxiliary data structure with working storage for RMHMC.
 */
typedef struct{
	double fx;			/* store current unormalised posterior density */
	double stepsize;	/* store current step size */
	double mL;		 	/* store num of leap-frog steps */
	int fIt;			/* store num of fixed point iterations */
	
	double* dfx;		/* store current likelihood gradient */
	double* cholMx;		/* store current cholesky factor of FI */
	double* invMx;		/* store current inverse of FI */
	double* dMx;		/* store current derivatives of FI */
	
	double target_acceptance;

	/* The rest is temporary working storage allocated in advance to 
	 * avoid repeated calls to malloc and free.
	 */
	double* new_x;		
	double  new_fx;     
	double* new_dfx;
	double* new_cholMx;
	double* new_invMx;
	double* new_dMx;
	
	
	double* momentum;
	double* new_momentum;
	double* tmpM;
	
	double* tr_invM_dM;
	double* atmp;
	double* btmp;
	double* p0;
	
} rmhmc_params;


/* implementation of the mcmc_kernel interface is in RMHMC_imp.c file 
   rmhmc_params type is used in RMHMC_imp.c so it must be defined before inclussion.
 */
#include "RMHMC_imp.c"


/* Free memory used by an rmhmc_params structure.
 * Arguments:
 *	params:		a pointer to the rmhmc_params structure.
 */
static void rmhmc_params_free(rmhmc_params* params){
	free(params->dfx);
	free(params->cholMx);
	free(params->invMx);
	free(params->dMx);
	free(params->new_x);
	free(params->new_dfx);
	free(params->new_cholMx);
	free(params->new_invMx);
	free(params->new_dMx);
	
	
	free(params->tr_invM_dM);
	free(params->tmpM);
	free(params->momentum);
	free(params->new_momentum);
	
	free(params->btmp);
	free(params->atmp);
	free(params->p0);
	
	free(params);
}

/* Allocate a new rmhmc_params structure for the working storage of RMHMC.
 * Arguments:
 *	N:				number of parameters, size of x.
 *	step_size:		initial step-size parameter.
 *	leap_frog:		initial number of leap frog steps.
 *	fixed_point:	number of fixed point iterations for the implicit solvers.
 * Returns:
 *	a new rmhmc_params structure.
 */
static rmhmc_params* rmhmc_params_alloc(int N, double step_size, double leap_frog, int fixed_point, double target_acceptance){
	
	rmhmc_params* params = (rmhmc_params*) malloc( sizeof(rmhmc_params) );
	if (params == 0) {
		/* TODO: write a proper error handler */
		fprintf(stderr,"rmhmc_params_alloc: malloc failed to allocate memory for rmhmc_params \n");
		return 0;		
	}
	/* initialise all pointers to 0 so we can safely call free in case of an error */
	params->dfx = 0; params->cholMx = 0; params->invMx = 0; params->dMx = 0; 
	
	params->new_x = 0; params->new_dfx = 0; params->new_cholMx = 0; 
	params->new_invMx = 0;  params->new_dMx = 0;
	
	params->momentum = 0; params->new_momentum = 0;
	params->tmpM = 0;

	params->tr_invM_dM = 0;
	params->atmp = 0; params->btmp = 0; params->p0 = 0;
	
	params->target_acceptance = target_acceptance; 

	/* =====  pre-allocate storage ======== */
	params->dfx = (double*) malloc( N * sizeof(double) );
	if (params->dfx == 0)
		goto error_cleanup;
		
	params->cholMx = (double*) malloc( N*N * sizeof(double) );
	if (params->cholMx == 0)
		goto error_cleanup;
	
	params->invMx = (double*) malloc( N*N * sizeof(double) );
	if (params->invMx == 0)
		goto error_cleanup;
	
	params->dMx = (double*) malloc( N*N*N * sizeof(double) );
	if (params->dMx == 0)
		goto error_cleanup;
	
	params->new_x = (double*) malloc( N * sizeof(double) );
	if (params->new_x == 0)
		goto error_cleanup;
	
	params->new_dfx = (double*) malloc( N * sizeof(double) );
	if (params->new_dfx == 0)
		goto error_cleanup;
		
	params->new_cholMx = (double*) malloc( N*N * sizeof(double) );
	if (params->new_cholMx == 0)
		goto error_cleanup;

	params->new_invMx = (double*) malloc( N*N * sizeof(double) );
	if (params->new_invMx == 0)
		goto error_cleanup;
	
	params->new_dMx = (double*) malloc( N*N*N * sizeof(double) );
	if (params->new_dMx == 0)
		goto error_cleanup;

	params->tr_invM_dM = (double*) malloc( N * sizeof(double) );
	if (params->tr_invM_dM == 0)
		goto error_cleanup;
	
	params->tmpM = (double*) malloc( N*N * sizeof(double) );
	if (params->tmpM == 0)
		goto error_cleanup;
	
	params->momentum = (double*) malloc( N * sizeof(double) );
	if (params->momentum == 0)
		goto error_cleanup;

	params->new_momentum = (double*) malloc( N * sizeof(double) );
	if (params->new_momentum == 0)
		goto error_cleanup;
	
	params->atmp = (double*) malloc( N * sizeof(double) );
	if (params->atmp == 0)
		goto error_cleanup;
	
	params->btmp = (double*) malloc( N * sizeof(double) );
	if (params->btmp == 0)
		goto error_cleanup;
	
	params->p0 = (double*) malloc( N * sizeof(double) );
	if (params->p0 == 0)
		goto error_cleanup;
	
		
	params->fx = 0;
	params->stepsize = step_size;
	params->mL = leap_frog;
	params->fIt = fixed_point;
	
	return(params);

error_cleanup:
	fprintf(stderr," rmhmc_params_alloc: malloc failed \n");
	rmhmc_params_free(params);
	return NULL;
}


/* Documentation in RMHMC.h */
rmhmc_model* rmhmc_model_alloc(fptrPosterior_rmhmc_all posteriorAll, 
							   fptrPosterior_rmhmc Metric,
							   fptrPrior_rnd Prx, void* model_params){
	
	rmhmc_model* model = (rmhmc_model*) malloc( sizeof (rmhmc_model) );
	if (model == 0){
		/* TODO: write a proper error handler */
		fprintf(stderr,"malloc failed to allocate memory for rmhmc_model\n");
		return 0;
	}
	
	model->PosteriorAll	= posteriorAll;
	model->Metric		= Metric;
	
	model->Prior_rnd = Prx;
	model->m_params = model_params;
	
	return model;
}

/* Documentation in RMHMC.h */
void rmhmc_model_free(rmhmc_model* model){
	free(model);
}

/* Documentation in RMHMC.h */
static void rmhmc_kernel_free(mcmc_kernel* kernel){
	rmhmc_params* params = (rmhmc_params*)kernel->kernel_params;
	rmhmc_params_free(params);
	free(kernel->x);
	gsl_rng* rng = (gsl_rng*)kernel->rng;
	gsl_rng_free(rng);
	free(kernel);
}

/* Documentation in RMHMC.h */
mcmc_kernel* rmhmc_kernel_alloc(int N, double step_size, double leap_frog, int fixed_point,
				rmhmc_model* model_function, unsigned long int seed, double target_acceptance){
	
	
  rmhmc_params* params = rmhmc_params_alloc(N, step_size, leap_frog, fixed_point, target_acceptance);
	if( params == 0 ){
		/* TODO: write a propper error handler */
		fprintf(stderr,"malloc failed to allocate memory for params in rmhmc_kernel_alloc \n");
		return 0;
	}
	
	mcmc_kernel* kernel = (mcmc_kernel*) malloc( sizeof(mcmc_kernel) );
	if (kernel == 0){
		/* TODO: write a proper error handler */
		rmhmc_params_free(params);
		fprintf(stderr,"malloc failed to allocate memory for mcmc_kernel in rmhmc_alloc.\n");
		return 0;
	}
	
	kernel->x = (double*) malloc( N * sizeof(double) );
	if (kernel->x == 0){
		/* TODO: write a proper error handler */
		rmhmc_params_free(params);
		free(kernel);
		fprintf(stderr,"malloc failed to allocate memory for mcnc_kernel.x in rmhmc_alloc \n");
		return 0;
	}
	
	gsl_rng_env_setup();
	gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
	if (r == 0){
		/* TODO: write a proper error handler */
		free(kernel->x);
		rmhmc_params_free(params);
		free(kernel);
		fprintf(stderr,"failed to create gsl_rng in rmhmc_alloc \n");
		return 0;
	}
	
	gsl_rng_set(r, seed);
	kernel->rng = r;
	
	kernel->N = N;
	kernel->model_function = model_function;
	kernel->kernel_params = params;
	
	/* rmhmc_kernel_* functions are defined in RMHMC_imp.c file */
	kernel->Sample = &rmhmc_kernel_sample;
	kernel->Adapt = &rmhmc_kernel_adapt;
	kernel->Init = &rmhmc_kernel_init;
	kernel->InitR = &rmhmc_kernel_init_rand;
	kernel->Free = &rmhmc_kernel_free;
	kernel->PrintStats = &rmhmc_kernel_print_stats;
	return kernel;
}
