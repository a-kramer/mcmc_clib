/*
 *  RMHMC_imp.c
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

/* Avoid double inclussion */

#ifndef __RMHMC_IMP_C__
#define __RMHMC_IMP_C__

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>

#include "RMHMC.h"

#define SWAP(a, b, tmp)  tmp = a; a = b; b = tmp;


/* Calculates trace terms used in Hamiltons equations.
 * Seccond term on the r.h.s. of equation (15) in Girolami and Calderhead (2011).
 * Arguments:
 *	N:		the number of parameters.
 *  invMx:	inverse of the metric ternsor evaluated at (x) in row major order.
 *  dMx:	partial derivatives of the metric tensor evaluated at (x) in row major order.
 * Result:
 *	tr_invM_dM: an array of N elements where 
 *				tr_invM_dM[i] = trace( M(x)^{-1} \frac{\partial M(x)}{\partial \theta_i ).
 */
static void calculateTraceTerms(int N, const double* invMx, const double* dMx, double* tr_invM_dM){
	
	int incIdx = N*N;
	int d,i,j;
	
	for (d = 0; d < N; d++) {
		/* initialise traces tr_invM_dM to 0 */
		tr_invM_dM[d] = 0.0;
		
		/* get Metric derivative w.r.t. for d^th parameter */
		const double* dMx_d = &dMx[d*incIdx];
		
		/* calulate trace terms using tr(A*B) = sum(sum( A.*B^T )) */
		for (i = 0; i < N; i++){
			for (j = 0; j < N; j++){
				tr_invM_dM[d] += invMx[i*N + j] *  dMx_d[j*N +i];
			}
		}
	}
	
}

/* Initialise RMHMC kernel with initial parameters x.
 * Arguments:
 *	kernel:		a pointer to the RMHMC kernel structure.
 *	x:			an array of N doubles with initial parameters. N must be 
 *				equal to kernel->N.
 * Result:
 * returns 0 for success and non-zero for error.
 */
static int rmhmc_kernel_init(mcmc_kernel* kernel, const double* x){
	int res,i,n;
	n = kernel->N;
	
	rmhmc_params* params = (rmhmc_params*)kernel->kernel_params;
	/* copy x to the kernel x state */
	for ( i=0; i < n; i++)
		kernel->x[i] = x[i];
	
	rmhmc_model* model = kernel->model_function;
	
	/* call user function to update all required quantities */
	res = model->PosteriorAll(x, model->m_params, &(params->fx), params->dfx, params->cholMx, params->dMx);
		
	/* TODO: write a proper error handler */
	if (res != 0){
		fprintf(stderr,"rmhmc_kernel_init: Likelihood function failed\n");
		return 1;
	}
	
	/* calculate cholesky factor for current metric */
	gsl_matrix_view cholMx_v = gsl_matrix_view_array(params->cholMx,n,n); 
	
	gsl_error_handler_t* old_handle =  gsl_set_error_handler_off();
	res = gsl_linalg_cholesky_decomp( &cholMx_v.matrix );
	if (res != 0){
		fprintf(stderr,"Error: matrix not positive definite in rmhmc_init.\n");
		return -1;
	}
	gsl_set_error_handler(old_handle);

	/* calculate inverse for current metric */
	gsl_matrix_view invMx_v = gsl_matrix_view_array(params->invMx,n,n);
	gsl_matrix_memcpy(&invMx_v.matrix, &cholMx_v.matrix );
	gsl_linalg_cholesky_invert(&invMx_v.matrix);
	
	/* calculate trace terms from equation (15) in Girolami and Calderhead (2011) */
	calculateTraceTerms(n, params->invMx, params->dMx, params->tr_invM_dM);
	
	return 0;
}

/* Initialise RMHMC kernel using random initial parameters sampled from the prior.
 * Arguments:
 *	kernel:		a pointer to the RMHMC kernel structure.
 * Result:
 *	returns 0 for success and non-zero for error.
 */
static int rmhmc_kernel_init_rand(mcmc_kernel* kernel){
	int n = kernel->N;
	
	gsl_rng* rng = (gsl_rng*) kernel->rng;
	rmhmc_model* model = kernel->model_function;
	
	/* sample random x from the prior */
	double x[n];		/* automatic alloc */
	model->Prior_rnd(rng, model->m_params, x);
	
	return rmhmc_kernel_init(kernel, x);
}

/* Adapt RMHMC parameters during the burn-in.
 * Arguments:
 *	kernel:		a pointer to the RMHMC kernel structure.
 *  acc_rate:   acceptance rate.
 */
static void rmhmc_kernel_adapt(mcmc_kernel* kernel, double acc_rate){
	rmhmc_params* params = (rmhmc_params*)kernel->kernel_params;

	/* Robins-Monro diminishing adaptation for the step-size */
	 params->stepsize  = params->stepsize - pow((double) kernel->t, - 0.2)*(params->target_acceptance - acc_rate);
	
	/*if (acc_rate > 0.8) {
		params->stepsize = 1.3*params->stepsize;
	}
	else if(acc_rate < 0.7){
		params->stepsize = 0.8*params->stepsize;
	}*/
	
	/* adapt number of leap-frog steps based on new step-size */
	double newL = 3.0/params->stepsize;
	if (newL > 30 )
		newL = 30;
	if (newL < 5)
		newL = 5;
	
	params->mL = newL;
}

/* Print statistics and information for the RMHMC sampler.
 * Arguments:
 *	kernel:		a pointer to the RMHMC kernel structure.
 *  outStream:  an output stream to write.
 */
static void rmhmc_kernel_print_stats(mcmc_kernel* kernel, FILE* outStream){
	rmhmc_params* params = (rmhmc_params*)kernel->kernel_params;
	fprintf(outStream, "Step-size: %5.4g LeapFrogSteps: %g \n",params->stepsize,params->mL);
	/* TODO: print additional sampler details */
}


/* Update momentum variables using an implicit solver for
 * equation (16) of Girolami and Calderhead (2011).
 * Arguments:
 *	state:		a pointer to internal working storage for RMHMC.
 *  p0:			current value of momentum parameters.
 *  iterations:	number of fixed point iterations. For iterations=1 then the function is 
 *				equivalent to the explicit solver for equation (18) of Girolami and Calderhead (2011).
 *	N:			number of parameters.
 *	stepSize:	integration step-size.
 *  Result:
 *	 The method directly updates the new_momentum array in the state structure.
 *	 returns number of iterations succussefuly completed.
 */
static int momentumNewtonUpdate(rmhmc_params* state, double* p0, int iterations, int N, double stepSize){
	
	gsl_vector_view new_p_v = gsl_vector_view_array(state->new_momentum, N);
	gsl_matrix_view new_invM_v = gsl_matrix_view_array(state->new_invMx, N, N);
	gsl_vector_view a_v = gsl_vector_view_array(state->atmp, N);
	gsl_vector_view b_v = gsl_vector_view_array(state->btmp, N);
	
	/* update trace terms */
	calculateTraceTerms(N, state->new_invMx, state->new_dMx, state->tr_invM_dM);
	
	int i,d;
	for (i = 0; i < iterations; i++) {
		
		/* a = invM*pNew; */
		gsl_blas_dgemv(CblasNoTrans, 1.0, &new_invM_v.matrix, &new_p_v.vector, 0.0, &a_v.vector);
		
		for (d = 0; d < N; d++) {
			gsl_matrix_view dMx_dv = gsl_matrix_view_array(&(state->new_dMx[d*N*N]), N,N);
			
			/* b =  dM/dx_n * invM*pNew  = dM/dx_n * a  */
			gsl_blas_dgemv (CblasNoTrans, 1.0, &dMx_dv.matrix, &a_v.vector, 0.0, &b_v.vector);
			
			/*p_update =  pNew^T*invM * dM/dx_n * invM*pNew = a^T * b  */
			double p_update = 0;
			gsl_blas_ddot(&a_v.vector, &b_v.vector, &p_update);
			
			p_update = p0[d] + 0.5*stepSize*( 0.5*p_update + state->new_dfx[d] - 0.5*state->tr_invM_dM[d] );
			state->new_momentum[d] = p_update;
		}
	}

	return i;
}

/* Update parameters using an implicit solver for
 * equation (17) of Girolami and Calderhead (2011).
 * Arguments:
 *	state:		a pointer to internal working storage for RMHMC.
 *  model:		a pointer to the rmhmc_model structure with pointers to user defined functions.
 *	N:			number of parameters.
 *	stepSize:	integration step-size.
 *  Result:
 *	 The method directly updates the new_x array in the state structure.
 *	 returns 0 for success or non-zero for failure.
 */
static int parametersNewtonUpdate(rmhmc_params* state, rmhmc_model* model, int N , double stepSize){
	
	gsl_vector_view new_x_v = gsl_vector_view_array(state->new_x, N);
	gsl_vector_view new_p_v = gsl_vector_view_array(state->new_momentum, N);
	gsl_matrix_view new_cholM_v = gsl_matrix_view_array(state->new_cholMx, N, N);
	
	/* temp copy of parameters */
	gsl_vector_view x0_v = gsl_vector_view_array(state->btmp, N);
	gsl_vector_memcpy(&x0_v.vector, &new_x_v.vector);
	
	/* temp copy of inverse Metric */
	gsl_matrix_view new_invM_v = gsl_matrix_view_array(state->new_invMx, N, N);
	gsl_matrix_view invM0_v = gsl_matrix_view_array(state->tmpM, N, N);
	gsl_matrix_memcpy(&invM0_v.matrix, &new_invM_v.matrix);
	
	gsl_vector_view a_v = gsl_vector_view_array(state->atmp, N);

	/* a = invM0*pNew */
	/* TODO: replace gsl_blas_dgemv with gsl_blas_dsymv since invM0_v.matrix is symetric */
	gsl_blas_dgemv(CblasNoTrans, 1.0, &invM0_v.matrix, &new_p_v.vector, 0.0, &a_v.vector);
	
	int iterations = state->fIt;
	int flag = 0;
	int i;
	for (i = 0; i < iterations; i++) {
		/* new_x = invM_new*p_new */
		/* TODO: replace gsl_blas_dgemv with gsl_blas_dsymv since inew_invM_v.matrix is symetric */
		gsl_blas_dgemv(CblasNoTrans, 1.0, &new_invM_v.matrix, &new_p_v.vector, 0.0, &new_x_v.vector);
		
		/* Calculates new_x_v = x0 + 0.5*stepSize*(invM_0*newP + newInvM*newP) */
		gsl_vector_add(&new_x_v.vector, &a_v.vector);
		gsl_vector_scale(&new_x_v.vector, 0.5*stepSize);
		gsl_vector_add(&new_x_v.vector, &x0_v.vector);
		
		/* calculate metric at the current position or update everything if this is the last iteration */
		if ( (i == iterations-1) )
			/* call user defined function for updating all quantities */
			model->PosteriorAll(state->new_x, model->m_params, &state->new_fx, state->new_dfx, state->new_cholMx, state->new_dMx);
		else
			/* call user defined function for updating only the metric ternsor */
			model->Metric(state->new_x, model->m_params, state->new_cholMx);
		
		/* calculate cholesky factor for current metric */
		gsl_error_handler_t* old_handle =  gsl_set_error_handler_off();
		flag = gsl_linalg_cholesky_decomp( &new_cholM_v.matrix );
		if (flag != 0){
			fprintf(stderr,"RMHMC: matrix not positive definite in parametersNewtonUpdate.\n");
			return flag;
		}
		gsl_set_error_handler(old_handle);
		
		/* calculate inverse for current metric */
		gsl_matrix_memcpy(&new_invM_v.matrix, &new_cholM_v.matrix );
		gsl_linalg_cholesky_invert(&new_invM_v.matrix);
	}
	return flag;
	
}

/* Calculate Hamiltonian energy equation (13) of Girolami and Calderhead (2011).
 * Arguments:
 *	N:			number of parameters.
 *  Lx:			unormalised posterior density at current parameter values.
 *	cholM:		cholesky factor of the metric tensor evaluated at current parameter values.
 *  invM:		inverse of metric tensor evaluated at current parameter values.
 *	p:			value of the momentum parameters.
 *  state:		a pointer to internal working storage for RMHMC.
 *  Result:
 *	 returns the log of the Hamiltonian energy.
 */
static double calculateHamiltonian(int N, double Lx, const double* cholM, const double* invM, 
								   const double* p, const rmhmc_params* state){
	double logDetM = 0;
	int d;
	
	gsl_vector_view a_v = gsl_vector_view_array(state->atmp, N);
	gsl_vector_const_view p_v = gsl_vector_const_view_array(p, N);
	gsl_matrix_const_view invM_v = gsl_matrix_const_view_array(invM, N, N);
	gsl_matrix_const_view cholM_v = gsl_matrix_const_view_array(cholM, N, N);

	//gsl_blas_dsymv(CblasUpper, 1.0, &invM_v.matrix, &p_v.vector, 0.0, &a_v.vector);
	gsl_blas_dgemv(CblasNoTrans, 1.0, &invM_v.matrix, &p_v.vector, 0.0, &a_v.vector);
	
	
	for (d = 0; d < N; d++)
		logDetM += log(gsl_matrix_get(&cholM_v.matrix, d, d));
	
	double quadTerm = 0;
	gsl_blas_ddot(&p_v.vector, &a_v.vector, &quadTerm);
	
	return Lx - logDetM - 0.5*quadTerm;
	
}

/* Copy variables from the previous sample in the RMHMC working storage to make 
 * space for the new proposal.
 * Arguments:
 *  state:		a pointer to internal working storage for RMHMC.
 *	kernel:	    a pointer to the RMHMC kernel data structure.
 *  Result:
 *	 returns zero for success and non-zero for failure.
 */
static int copyStateVariables(rmhmc_params* state, mcmc_kernel* kernel){
	int N = kernel->N;
	gsl_vector_view p_v = gsl_vector_view_array(state->momentum, N);
	gsl_vector_view new_p_v = gsl_vector_view_array(state->new_momentum, N);
	gsl_vector_memcpy(&new_p_v.vector, &p_v.vector);
	
	gsl_vector_view x_v = gsl_vector_view_array(kernel->x, N);
	gsl_vector_view new_x_v = gsl_vector_view_array(state->new_x, N);
	gsl_vector_memcpy(&new_x_v.vector, &x_v.vector);
	
	gsl_vector_view dfx_v = gsl_vector_view_array(state->dfx, N);
	gsl_vector_view new_dfx_v = gsl_vector_view_array(state->new_dfx, N);
	gsl_vector_memcpy(&new_dfx_v.vector, &dfx_v.vector);
		
	gsl_matrix_view invM_v = gsl_matrix_view_array(state->invMx, N, N);
	gsl_matrix_view new_invM_v = gsl_matrix_view_array(state->new_invMx, N, N);
	gsl_matrix_memcpy(&new_invM_v.matrix, &invM_v.matrix);
	
	int i;
	for(i = 0; i < N*N*N; i++){ 
		state->new_dMx[i]  = state->dMx[i];
	}
	
	return 0;
}

/* Draw one sample from the Markov Chain using the RMHMC algorithm.
 * Arguments:
 *	kernel:	    a pointer to the RMHMC kernel data structure.
 *  Result:
 *	 returns zero for success and non-zero for failure.
 *	 the new sample is directly updated in kernel->x.
 *	 acc is set to 0 if the chain  in the previous state (reject)
 *   and 1 if the chain made a transition to a new state (accept)
 */
static int rmhmc_kernel_sample(mcmc_kernel* kernel, int* acc){

	rmhmc_params* state = (rmhmc_params*) kernel->kernel_params;
	rmhmc_model* model	= (rmhmc_model*) kernel->model_function;
	gsl_rng* rng = (gsl_rng*) kernel->rng;
	
	int N = kernel->N;
	double stepSize = state->stepsize;
	int fIt = state->fIt;
	
	/* sample momentum variables from multivariate normal with covariance Mx */
	double* p = state->momentum;
	int d;
	for (d = 0; d < N; d++)
		p[d] = gsl_ran_ugaussian(rng);
	
	gsl_vector_view p_v = gsl_vector_view_array(state->momentum, N);
	gsl_matrix_view cholM_v = gsl_matrix_view_array(state->cholMx, N, N);
	/* p = cholM*p */
	gsl_blas_dtrmv (CblasUpper, CblasTrans, CblasNonUnit, &cholM_v.matrix, &p_v.vector);
	
	/* randomise direction of integration */
	double randDir = gsl_rng_uniform(rng);
	if (randDir > 0.5) 
		stepSize = -1.0*stepSize;
	
	/* randomise number of leap-frog steps */
	int L = 1 + gsl_rng_uniform_int(rng, state->mL);
	
	/* Generalised leap-frog integrator */
	copyStateVariables(state, kernel);

	gsl_vector_view new_p_v = gsl_vector_view_array(state->new_momentum, N);
	gsl_vector_view tmpVect = gsl_vector_view_array(state->p0, N);
	
	int l;
	int flag = 0;
	for (l = 0; l < L; l++) {
		
		/* momentum Newton update */
		/* temp copy of momentum variables */
		gsl_vector_memcpy(&tmpVect.vector, &new_p_v.vector);
		momentumNewtonUpdate(state, state->p0, fIt, N, stepSize);
		
		/* parameters Newton update */
		flag = parametersNewtonUpdate(state, model, N, stepSize);
		
		if (flag != 0){
			fprintf(stderr,"RMHMC: Error in parameter Newton update. Reject step.\n");
			*acc = 0;
			return flag;
		}
		
		/* single Newton update step for momentum variables */
		momentumNewtonUpdate(state, state->new_momentum, 1, N, stepSize);
		
	}
	
	/* calculate Hamiltonian energy for current state */
	double H_c = calculateHamiltonian(N, state->fx, state->cholMx, state->invMx, state->momentum, state);
	
	/* calculate Hamiltonian energy for proposed state */
	double H_p = calculateHamiltonian(N, state->new_fx, state->new_cholMx, state->new_invMx, state->new_momentum, state);
	
	/* Accept/reject using Metropolis-Hastings ratio */
	double mh_ratio = H_p - H_c;
	double rand_dec = log(gsl_rng_uniform(rng));
	
	if ( (mh_ratio > 0.0)||(mh_ratio > rand_dec) ) {
		*acc = 1;
		double* tmp;
		SWAP(kernel->x, state->new_x, tmp);
		SWAP(state->dfx, state->new_dfx, tmp);
		SWAP(state->cholMx, state->new_cholMx, tmp);
		SWAP(state->invMx, state->new_invMx, tmp);
		SWAP(state->dMx, state->new_dMx, tmp);
		state->fx  = state->new_fx;
		
	}else {
		
		*acc = 0;
	}
	
	return 0;
	
}




#endif
