/*
 *  mv_norm.h
 *  Incomplete interface for multivariate normal densities
 *
 */

#ifndef __MVNORM_H__
#define __MVNORM_H__

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>	

/* multivariate normal log density function    */
/*
 *    x       random vetor of size n
 *    mean    vector of means of size n. Pass NULL if mean is zero.
 *    cholPr  upper triangular cholesky of the precision matrix of dimension n x n
*/
double log_mv_norm_pdf_cholP(const gsl_vector* x, const gsl_vector* m, const gsl_matrix* cholPr);


/* multivariate normal random variates         */
/*
 *	  rng     random number generator
 *    mean    vector of means of size n, Pass NULL if mean is zero.
 *    cholPr  upper triangular cholesky of the precision matrix of dimension n x n
 *	  result  vector of random variates of size n
 */
	
int mnv_norm_rnd_cholPr(gsl_rng* rng, const gsl_vector* m, const gsl_matrix* cholPr,
							gsl_vector* result);

/* multivariate normal log density function    */
/*
 *    x       random vetor of size n
 *    mean    vector of means of size n, Pass NULL if mean is zero.
 *    var     covariance matrix of dimension n x n
*/
double log_mv_norm_pdf_var(const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var);

	
/* multivariate normal log density function    */
/*
 *    x       random vetor of size n
 *    mean    vector of means of size n, Pass NULL if mean is zero.
 *    cholVr  lower triangular cholesky of the covariance matrix of dimension n x n
*/
	
double log_mv_norm_pdf_cholV(const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *cholVr);


#ifdef __cplusplus
}
#endif

#endif