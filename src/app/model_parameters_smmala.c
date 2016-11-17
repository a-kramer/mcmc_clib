#include "model_parameters_smmala.h"

int ode_model_parameters_alloc(ode_model_parameters *omp, const problem_size *ps){
  /* problem_size structure contains:
   * N: number of states
   * D: number of model parameters
   * F: number of output functions
   * U: number of input parameters: u[j*U+i]: i=0,...,U-1
   * C: number of different experimental conditions (No. of inputs u): u[j*U+i]: j=0,...,C-1
   * T: number of measured time points
   * n1,n2 the size of the normalisation array 
   * (time and state indices for normalisation):
   *  e.g., for n1=2; n2=F=3,
   *  [normalisation]
   *   2,1,2
   *   0,3,1
   *  [/normalisation]
   * means that output function 0 is normalised by 
   * output function 2 at its initial condition.
   * Output function 1 is normalised by itself at time point t[3]
   * Output function 2 is also normalised by itself but at time point t[1]
   */
  
  gsl_vector **y,**fy;
  gsl_matrix **yS,**fyS;

  int N=ps->N;
  int C=ps->C;
  int F=ps->F;
  int D=ps->D;
  int U=ps->U;
  int T=ps->T;
  int P=ps->P;
  
  printf("# allocating model parameters with: N=%i,\tP=%i,\tT=%i.\n",N,P,T);

  /* initialise gsl_vectors and matrices these will hold the
   * differential equation data for each experimental condition c and
   * time point t_j, like this: y[c*T+j]=gsl_vector(N)
   */ 

  omp->normalisation=gsl_matrix_alloc(ps.n1,ps.n2);
  
  omp->y=(gsl_vector*) malloc(sizeof(gsl_vector*)*C*T);
  omp->fy=(gsl_vector*) malloc(sizeof(gsl_vector*)*C*T);
  omp->yS=(gsl_matrix*) malloc(sizeof(gsl_matrix*)*C*T);
  omp->fyS=(gsl_matrix*) malloc(sizeof(gsl_matrix*)*C*T);
  omp->oS=(gsl_matrix*) malloc(sizeof(gsl_matrix*)*C*T);
    y=omp->y;
   fy=omp->fy;
   yS=omp->yS;
  fyS=omp->fyS;

  for (i=0;i<C*T;i++){
    y[i]=gsl_vector_alloc(N);
   fy[i]=gsl_vector_alloc(F);
   yS[i]=gsl_matrix_alloc(P,N);
  fyS[i]=gsl_matrix_alloc(P,F);
   oS[i]=gsl_matrix_alloc(D,F);
  }
  /* during burn-in, we slowly increase beta from 0 to 1; if no
   * burn-in is performed, beta needs to be 1.0
   */
  omp->beta=1.0;
  
  omp->reference_y=gsl_matrix_alloc(T,N);
  omp->output_C=gsl_matrix_alloc(F,N);

  // [\rho=exp(theta), u]
  omp->exp_x_u=gsl_vector_alloc(P);
  // time
  omp->t=gsl_vector_alloc(T);

  // Data (t is the fast index, c the slow one). One huge block of Data:
  omp->Data=gsl_matrix_alloc(C*T,F);
  omp->sdData=gsl_matrix_alloc(C*T,F);
  // sometimes it is convenient to do vector operations on single
  // lines of the data.
  // we create some row views for that purpose:
  omp->data_row=(gsl_vector_view*) malloc(sizeof(gsl_vector_view)*C*T);
  omp->sd_data_row=(gsl_vector_view*) malloc(sizeof(gsl_vector_view)*C*T);
  // list of vector shorcuts to the views, no need to allocate the
  // vectors themselves
  omp->data=(gsl_vector**) malloc(sizeof(gsl_vector*)*C*T);
  omp->sd_data=(gsl_vector**) malloc(sizeof(gsl_vector*)*C*T);
  for (i=0;i<C*T;i++){
    omp->data_row[i]=gsl_matrix_row(omp->Data,i);
    omp->data[i]=&(omp->data_row[i]);
    omp->sd_data_row[i]=gsl_matrix_row(omp->Data,i);
    omp->sd_data[i]=&(omp->data_row[i]);
  }
  

  
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
  //omp->yS=gsl_matrix_alloc(P,N);
  omp->reference_yS=gsl_matrix_alloc(T*P,N);
  //omp->fyS=gsl_matrix_alloc(P,F);
  omp->reference_fyS=gsl_matrix_alloc(T*P,F);
  //omp->oS=gsl_matrix_alloc(P,F);

  return EXIT_SUCCESS;
}

int ode_model_parameters_free(ode_model_parameters *omp){
  int C=omp->input_u->size1;
  int T=omp->t->size;
  int CT=C*T;
  int i;
  for (i=0;i<CT;i++){
    gsl_vector_free(omp->y[i]);
    gsl_vector_free(omp->fy[i]);
    gsl_matrix_free(omp->yS[i]);
    gsl_matrix_free(omp->fyS[i]);
    gsl_matrix_free(omp->oS[i]);
  }
  free(omp->data_row);
  free(omp->sd_data_row);
  free(omp->data);
  free(omp->data);
  
  gsl_matrix_free(omp->normalisation);
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
  //gsl_matrix_free(omp->oS);

  return EXIT_SUCCESS;
 
}
