#include "model_parameters_rmhmc.h"

//ode_model_parameters_alloc(omp, D, N, F, T, U, C);
int ode_model_parameters_alloc(ode_model_parameters *omp, int D, int N, int F,  int T, int U, int C){
  int cvode_N=omp->solver->odeModel->N;
  int P=D+U;
  printf("allocating model parameters with: N=%i,\tP=%i,\tT=%i.\n",N,P,T);
  
  omp->cvodes_state = (double*) malloc(sizeof(double)*cvode_N);
  omp->cvodes_reference_state = (double*) malloc(sizeof(double)*T*cvode_N);

  omp->fy=gsl_vector_alloc(F);
  omp->fyS=gsl_matrix_alloc(P,F);
  omp->dfyS=gsl_matrix_alloc(P*P,F);

  omp->rfy=gsl_vector_alloc(F);
  omp->rfyS=gsl_matrix_alloc(P,F);
  omp->drfyS=gsl_matrix_alloc(P*P,F);

  omp->oS=gsl_matrix_alloc(D,F);
  omp->doS=gsl_matrix_alloc(D*D,F);

  omp->output_C=gsl_matrix_alloc(F,N);
  omp->exp_x_u=gsl_vector_alloc(P);

  // time
  omp->t=gsl_vector_alloc(T);

  // Data
  omp->Data=gsl_matrix_alloc(T*C,F);
  omp->sdData=gsl_matrix_alloc(T*C,F);

  // inputs
  omp->input_u=gsl_matrix_alloc(C,U);
  omp->reference_u=gsl_vector_alloc(U);

  // prior 
  omp->prior_tmp_a=gsl_vector_alloc(D);
  omp->prior_tmp_b=gsl_vector_alloc(D);
  omp->prior_mu=gsl_vector_alloc(D);
  omp->prior_inverse_cov=gsl_matrix_alloc(D,D);

  omp->oS=gsl_matrix_alloc(P,F);

  return EXIT_SUCCESS;
}

int ode_model_parameters_free(ode_model_parameters *omp){

  free(omp->cvodes_state);
  free(omp->cvodes_reference_state);

  gsl_matrix_free(omp->output_C);
  gsl_vector_free(omp->fy);
  gsl_matrix_free(omp->fyS);
  gsl_matrix_free(omp->dfyS);

  gsl_vector_free(omp->rfy);
  gsl_matrix_free(omp->rfyS);
  gsl_matrix_free(omp->drfyS);

  gsl_matrix_free(omp->oS);
  gsl_matrix_free(omp->doS);

  gsl_vector_free(omp->exp_x_u);
  // time
  gsl_vector_free(omp->t);

  // Data
  gsl_matrix_free(omp->Data);
  gsl_matrix_free(omp->sdData);

  // inputs
  gsl_matrix_free(omp->input_u);
  gsl_vector_free(omp->reference_u);

  // prior 
  gsl_vector_free(omp->prior_tmp_a);
  gsl_vector_free(omp->prior_tmp_b);
  gsl_vector_free(omp->prior_mu);
  gsl_matrix_free(omp->prior_inverse_cov);

  // sensitivities
  gsl_matrix_free(omp->oS);

  return EXIT_SUCCESS;
 
}
