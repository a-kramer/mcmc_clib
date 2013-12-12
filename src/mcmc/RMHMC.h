/*
 *  RMHMC.h
 *
 *  Copyright (C) 2012 Vassilios Stathopoulos stathv@gmail.com
 *
 *	This program is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	Foobar is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef __RMHMC_H__
#define __RMHMC_H__

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
	
#include <gsl/gsl_rng.h>
#include "mcmc_kernel.h"
	
	/* User defined function for caclulating the metric ternsor.
	 * Arguments:
	 *	x:				values for the model parameters.
	 *  model_params:   a pointer to a user defined structure with aditional model data and working storage.
	 * Returns:
	 *  FI:				a double aray with size N*N with elements of the Metric tensor stored in row major order.
	 *					N the size of x.
	 */
	typedef int (*fptrPosterior_rmhmc)(const double* x,  void* model_params, double* FI);
	
	/* User defined function for caclulating the unormalised posterior, gradient, metric ternsor and partial
	 * derivatives of the metric tensor.
	 * Arguments:
	 *	x:				values for the model parameters.
	 *  model_params:   a pointer to a user defined structure with aditional model data and working storage.
	 * Returns:
	 *  fx:				unormalised posterior density, p(DATA|x)p(x)
	 *  dfx:			gradient of the unormalised posterior density.
	 *  FI:				a double aray with size N*N with elements of the Metric tensor stored in row major order.
	 *					N the size of x.
	 *	dFI:            a double array with size N*N*N with partial derivatives of the Metric tensor w.r.t. x.
	 *					The i,j element of the k^th partial derivative matrix is stored at dFI[k*N*N + i*N + j].
	 */
	typedef int (*fptrPosterior_rmhmc_all)(const double* x,  void* model_params, double* fx, 
										   double* dfx, double* FI, double* dFI );

	/* Auxiliary data structure to hold pointers to user defined functions and model axuliary working storage and 
	 * aditional parameters or constants.
	 */
	typedef struct{
		fptrPosterior_rmhmc_all PosteriorAll;
		fptrPosterior_rmhmc Metric;		
		fptrPrior_rnd  Prior_rnd;
		void* m_params;
	} rmhmc_model;
	
	
	/* Allocates a new model structure for the RMHMC sampler.
	 * Arguments:
	 *   posteriorAll:	user defined function for calculating all nessesary quantities.
	 *	 Metric:		user defined function for calculating the metric tensor.
	 *   Prx:			user defined function for drawing random samples from the prior.
	 *	 model_params:  user defined structure for aditional model parameters and working storage.
	 * Returns:
	 * a new rmhmc_model data structure.
	 */
	rmhmc_model* rmhmc_model_alloc(fptrPosterior_rmhmc_all posteriorAll, 
								   fptrPosterior_rmhmc Metric,
								   fptrPrior_rnd Prx, void* model_params);
	
	/* De-allocates an rmhmc_model structure.
	 * Arguments:
	 *   model:		a pointer to the rmhmc_model structure to dealocate.
	 */
	void rmhmc_model_free(rmhmc_model* model);

	
	/* Allocates a new MCMC kernel structure for the RMHMC sampler.
	 * Arguments:
	 *   N:					number of parameters, size of x.
	 *	 step_size:			initial step-size parameter.
	 *   leap_frog:			initial number of leap-frog steps.
	 *   fixed_point:		number of fixed point iterations for implicit solvers.
	 *	 model_function:	a pointer to the data structure with the user defined functions.
	 * Returns:
	 * a new mcmc_kernel data structure for RMHMC.
	 */	
	mcmc_kernel* rmhmc_kernel_alloc(int N, double step_size, double leap_frog, int fixed_point,
					rmhmc_model* model_function, unsigned long int seed, double target_acceptance);
		
		
#ifdef __cplusplus
	}
#endif
	
#endif
