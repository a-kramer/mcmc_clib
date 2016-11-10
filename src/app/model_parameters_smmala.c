#include "model_parameters_smmala.h"

int ode_model_parameters_alloc(ode_model_parameters *omp, int D, int N, int F,  int T, int U, int C){
  /* N: number of states
   * D: number of model parameters
   * F: number of output functions
   * U: number of input parameters: u[j*U+i]: i=0,...,U-1
   * C: number of different experimental conditions (No. of inputs u): u[j*U+i]: j=0,...,C-1
   * T: number of measured time points
   */

  /* the Data structure: concatenation
   * row index: 0,...,P-1
   * column index: 0,...,T*U-1
   * D=cat(2,Data(u_1),Data(u_2),...,Data(u_U))
   */
  
  gsl_vector **y,**fy;
  gsl_matrix **yS,**fyS;

  int P=omp->solver->odeModel->P;
  printf("# allocating model parameters with: N=%i,\tP=%i,\tT=%i.\n",N,P,T);

  /* initialise gsl_vectors and matrices these will hold the
   * differential equation data for each experimental condition c and
   * time point t_j, like this: y[c*T+j]=gsl_vector(N)
   */ 
  
  omp->y=(gsl_vector*) malloc(sizeof(gsl_vector*)*C*T);
  omp->fy=(gsl_vector*) malloc(sizeof(gsl_vector*)*C*T);
  omp->yS=(gsl_matrix*) malloc(sizeof(gsl_matrix*)*C*T);
  omp->fyS=(gsl_matrix*) malloc(sizeof(gsl_matrix*)*C*T);

    y=omp->y;
   fy=omp->fy;
   yS=omp->yS;
  fyS=omp->fyS;

  for (i=0;i<C*T;i++)   y[i]=gsl_vector_alloc(N);
  for (i=0;i<C*T;i++)  fy[i]=gsl_vector_alloc(N);
  for (i=0;i<C*T;i++)  yS[i]=gsl_matrix_alloc(D,N);
  for (i=0;i<C*T;i++) fyS[i]=gsl_matrix_alloc(D,N);

  
  omp->beta=1.0;

  
  gsl_vector_alloc(N);
  omp->reference_y=gsl_matrix_alloc(T,N);
  omp->output_C=gsl_matrix_alloc(F,N);

  // [\rho=exp(theta), u]
  omp->exp_x_u=gsl_vector_alloc(P);
  // time
  omp->t=gsl_vector_alloc(T);

  // Data
  omp->Data=gsl_matrix_alloc(T*C,F);
  omp->sdData=gsl_matrix_alloc(T*C,F);

  // outputs
  omp->fy=gsl_vector_alloc(F);
  omp->reference_fy=gsl_matrix_alloc(T,F);

  // inputs
  omp->input_u=gsl_matrix_alloc(C,U);
  omp->reference_u=gsl_vector_alloc(U);

  // prior 
  omp->prior_tmp_a=gsl_vector_alloc(D);
  omp->prior_tmp_b=gsl_vector_alloc(D);
  omp->prior_mu=gsl_vector_alloc(D);
  omp->prior_inverse_cov=gsl_matrix_alloc(D,D);

  // sensitivities

  omp->yS0=gsl_matrix_calloc(P,N); //yS0=zeros(P,N);
  omp->yS=gsl_matrix_alloc(P,N);
  omp->reference_yS=gsl_matrix_alloc(T*P,N);
  omp->fyS=gsl_matrix_alloc(P,F);
  omp->reference_fyS=gsl_matrix_alloc(T*P,F);
  omp->oS=gsl_matrix_alloc(P,F);

  return EXIT_SUCCESS;
}

int ode_model_parameters_free(ode_model_parameters *omp){
  int CT=omp->Data->size1;
  int i;
  for (i=0;i<CT;i++) gsl_vector_free(omp->y[i]);
  for (i=0;i<CT;i++) gsl_vector_free(omp->fy[i]);
  for (i=0;i<CT;i++) gsl_matrix_free(omp->yS[i]);
  for (i=0;i<CT;i++) gsl_matrix_free(omp->fyS[i]);


  gsl_matrix_free(omp->reference_y);
  gsl_matrix_free(omp->output_C);
  gsl_vector_free(omp->exp_x_u);
  
  // time
  gsl_vector_free(omp->t);

  // Data
  gsl_matrix_free(omp->Data);
  gsl_matrix_free(omp->sdData);

  // outputs
  gsl_matrix_free(omp->reference_fy);

  // inputs
  gsl_matrix_free(omp->input_u);
  gsl_vector_free(omp->reference_u);

  // prior 
  gsl_vector_free(omp->prior_tmp_a);
  gsl_vector_free(omp->prior_tmp_b);
  gsl_vector_free(omp->prior_mu);
  gsl_matrix_free(omp->prior_inverse_cov);

  // sensitivities
  gsl_matrix_free(omp->yS0);
  gsl_matrix_free(omp->reference_yS);
  gsl_matrix_free(omp->reference_fyS);
  gsl_matrix_free(omp->oS);

  return EXIT_SUCCESS;
 
}
