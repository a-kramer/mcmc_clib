/*
 *  mv_norm.c
 *  
 *
 *  Created by Vassilios Stathopoulos on 15/09/2011.
 *  Copyright 2011 Computing Science. All rights reserved.
 *
 */

#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include "mv_norm.h"


static const double log2Pi = 1.837877066409345;

double log_mv_norm_pdf_cholP(const gsl_vector* x, const gsl_vector* m, const gsl_matrix* cholPr){

	
	int i,n;
	double qdr,logDet,lpdf;

	n = x->size;
	
	double work[n]; /* automatic allocation is faster than dynamic */
	gsl_vector_view work_v = gsl_vector_view_array(work, n);
	gsl_vector* work_mv = &work_v.vector;
	
	/* p = sum(log(diag(L^{T}))) -0.5 sum( (L^{T}*(x-m)).^2 )/stepsize 
	 where L the lower triangular chol of H */	
	gsl_vector_memcpy(work_mv, x);
	if (m != NULL)
		gsl_vector_sub(work_mv, m);
	
	gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, cholPr, work_mv);
	gsl_blas_ddot(work_mv,work_mv, &qdr);
	
	logDet = 0.0;
	for(i=0; i<n; i++)
		logDet += log(gsl_matrix_get(cholPr,i,i));
	
	lpdf = - 0.5*n*log2Pi + logDet - 0.5*qdr;
	
	
	return lpdf;
}

int mnv_norm_rnd_cholPr(gsl_rng* rng, const gsl_vector* m, const gsl_matrix* cholPr,
						gsl_vector* result){
	/* TODO: error checking and return error code */

	int i,n;
	n = cholPr->size1;
	
	for(i=0; i<n; i++)
		gsl_vector_set( result, i, gsl_ran_ugaussian(rng) );
	
	/* For covariance A then L*randn, for covariance A^{-1} then L^{-T}*rand */
	gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, cholPr, result);
	/*gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, cholPr, result);*/
	if (m != NULL)
		gsl_vector_add(result, m);
	
	return 0;
}


double log_mv_norm_pdf_var(const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var){
	
	int n = x->size;

	double chol_arr[n*n]; /* automatic allocation is faster than dynamic */
	double x0_arr[n];
	gsl_matrix_view chol_v	= gsl_matrix_view_array(chol_arr, n, n);
	gsl_matrix* chol		= &chol_v.matrix;
	gsl_vector_view x0_v	= gsl_vector_view_array(x0_arr, n);
	gsl_vector* x0			= &x0_v.vector;

	gsl_vector_memcpy(x0, x);
	if (mean != NULL)
		gsl_vector_sub(x0, mean);					/* x0 = x - m */
	
	gsl_matrix_memcpy( chol, var );
	gsl_linalg_cholesky_decomp( chol );
	gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, chol, x0);
	/*gsl_linalg_cholesky_svx(chol, x0);		// x0 = (x - m) / V */
	double qdr;
	gsl_blas_ddot(x0, x0, &qdr);			/* qdr = (x-m)^T V^-1 (x-m) */
	
	double logDet = 0.0;
	int i;
	for (i = 0; i < n; i++)
		logDet += log(gsl_matrix_get(chol,i,i));

	double log_pdf = - 0.5*n*log2Pi - logDet - 0.5*qdr;
	return log_pdf;
}


double log_mv_norm_pdf_cholV(const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *cholVr){
	
	int n = x->size;
	
	double x0_arr[n];	/* automatic allocation is faster than dynamic */
	gsl_vector_view x0_v	= gsl_vector_view_array(x0_arr, n);
	gsl_vector* x0			= &x0_v.vector;
	
	gsl_vector_memcpy(x0, x);
	if(mean != NULL )
		gsl_vector_sub(x0, mean);					/* x0 = x - m */
	
	gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, cholVr, x0);
	/*gsl_linalg_cholesky_svx(chol, x0);		// x0 = (x - m) / V */
	double qdr;
	gsl_blas_ddot(x0, x0, &qdr);			/* qdr = (x-m)^T V^-1 (x-m) */
	
	double logDet = 0.0;
	int i;
	for (i = 0; i < n; i++)
		logDet += log(gsl_matrix_get(cholVr,i,i));
	
	double log_pdf = - 0.5*n*log2Pi - logDet - 0.5*qdr;
	return log_pdf;
}
