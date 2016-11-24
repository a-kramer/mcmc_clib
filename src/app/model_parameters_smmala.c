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
   * this can repeat for every experimental condition.
   */
  gsl_vector **r_y,**r_fy;
  gsl_matrix **r_yS, **r_fyS;
  gsl_vector **y,**fy;
  gsl_matrix **yS,**fyS, **oS;  
  int i;
  int N=ps->N; // number of state variables
  int C=ps->C; // number of experimental conditions
  int F=ps->F; // number of output functions
  int D=ps->D; // number of sampling parameters (to be estimated)
  int U=ps->U; // number of input parameters (known)
  int T=ps->T; // number of measurement time instances
  int P=ps->P; // number of total parameters D+U (a consistency check
	       // between ode_model and mcmc configuration file)
  
  //printf("# allocating model parameters with: N=%i,\tP=%i,\tT=%i.\n",N,P,T);
  /* initialise gsl_vectors and matrices these will hold the
   * differential equation data for each experimental condition c and
   * time point t_j, like this: y[c*T+j]=gsl_vector(N)
   */
  omp->tmpF=gsl_vector_alloc(F);
  //printf("alloc norm f and t\n");
  omp->norm_f=gsl_matrix_int_alloc(C,F);
  omp->norm_t=gsl_matrix_int_alloc(C,F);
  //printf("array alloc y, fy, yS, fyS, oS\n");
  
  omp->y=(gsl_vector**) malloc(sizeof(gsl_vector*)*C*T);
  omp->fy=(gsl_vector**) malloc(sizeof(gsl_vector*)*C*T);
  omp->yS=(gsl_matrix**) malloc(sizeof(gsl_matrix*)*C*T);
  omp->fyS=(gsl_matrix**) malloc(sizeof(gsl_matrix*)*C*T);
  omp->oS=(gsl_matrix**) malloc(sizeof(gsl_matrix*)*C*T);
  
    y=omp->y;
   fy=omp->fy;
   yS=omp->yS;
  fyS=omp->fyS;
   oS=omp->oS;
   //printf("gsl alloc y, fy, yS, fyS, oS\n");

  for (i=0;i<C*T;i++){
    y[i]=gsl_vector_alloc(N);
   fy[i]=gsl_vector_alloc(F);
   yS[i]=gsl_matrix_alloc(P,N);
  fyS[i]=gsl_matrix_alloc(P,F);
   oS[i]=gsl_matrix_alloc(D,F);
  }
  //printf("ptr alloc reference: y, fy, yS, fyS\n");

  omp->reference_y=(gsl_vector**) malloc(sizeof(gsl_vector*)*T);
  omp->reference_fy=(gsl_vector**) malloc(sizeof(gsl_vector*)*T);
  omp->reference_yS=(gsl_matrix**) malloc(sizeof(gsl_matrix*)*T);
  omp->reference_fyS=(gsl_matrix**) malloc(sizeof(gsl_matrix*)*T);

  r_y=omp->reference_y;
  r_fy=omp->reference_fy;
  r_yS=omp->reference_yS;
  r_fyS=omp->reference_fyS;
  //printf("gsl alloc reference: y, fy, yS, fyS\n");

  for (i=0;i<T;i++){
    r_y[i]=gsl_vector_alloc(N);
    r_fy[i]=gsl_vector_alloc(F);
    r_yS[i]=gsl_matrix_alloc(P,N);
    r_fyS[i]=gsl_matrix_alloc(P,F);
  }
  /* during burn-in, we slowly increase beta from 0 to 1; if no
   * burn-in is performed, beta needs to be 1.0
   */
  omp->beta=1.0;
  
  // [\rho=exp(theta), u]
  //printf("gsl alloc exp_x_u\n");

  omp->exp_x_u=gsl_vector_alloc(P);
  // time
  printf("gsl alloc t\n");
  omp->t=gsl_vector_alloc(T);

  // Data (t is the fast index, c the slow one). One huge block of Data:
  //printf("gsl alloc Data\n");

  omp->Data=gsl_matrix_alloc(C*T,F);
  omp->sdData=gsl_matrix_alloc(C*T,F);
  // sometimes it is convenient to do vector operations on single
  // lines of the data.
  // we create some row views for that purpose:
  //printf("ptr alloc Data row views\n");

  omp->data_row=(gsl_vector_view*) malloc(sizeof(gsl_vector_view)*C*T);
  omp->sd_data_row=(gsl_vector_view*) malloc(sizeof(gsl_vector_view)*C*T);
  // list of vector shorcuts to the views, no need to allocate the
  // vectors themselves
  omp->data=(gsl_vector**) malloc(sizeof(gsl_vector*)*C*T);
  omp->sd_data=(gsl_vector**) malloc(sizeof(gsl_vector*)*C*T);
  //printf("assign Data row views\n");

  for (i=0;i<C*T;i++){
    omp->data_row[i]=gsl_matrix_row(omp->Data,i);
    omp->data[i]=&(omp->data_row[i].vector);
    omp->sd_data_row[i]=gsl_matrix_row(omp->sdData,i);
    omp->sd_data[i]=&(omp->sd_data_row[i].vector);
  }
  // shortcuts for input rows;
  //printf("ptr alloc input row views\n");

  // inputs
  omp->input_u=gsl_matrix_alloc(C,U);
  omp->reference_u=gsl_vector_alloc(U);
  omp->input_u_row=(gsl_vector_view*) malloc(sizeof(gsl_vector_view)*C);
  omp->u=(gsl_vector **) malloc(sizeof(gsl_vector*)*C);
  //printf("assign input row views\n");
  for (i=0;i<C;i++){
    omp->input_u_row[i]=gsl_matrix_row(omp->input_u,i);
    omp->u[i]=&(omp->input_u_row[i].vector);
  }
  

  // prior 
  omp->prior_tmp_a=gsl_vector_alloc(D);
  omp->prior_tmp_b=gsl_vector_alloc(D);
  omp->prior_mu=gsl_vector_alloc(D);
  omp->prior_inverse_cov=gsl_matrix_alloc(D,D);

  // sensitivities

  omp->yS0=gsl_matrix_calloc(P,N); //yS0=zeros(P,N);
  return EXIT_SUCCESS;
}

int ode_model_parameters_free(ode_model_parameters *omp){
  int C=omp->input_u->size1;
  int T=omp->t->size;
  int CT=C*T;
  int i;


  gsl_matrix_free(omp->yS0);
  for (i=0;i<T;i++){
    gsl_vector_free(omp->reference_y[i]);
    gsl_vector_free(omp->reference_fy[i]);
    gsl_matrix_free(omp->reference_yS[i]);
    gsl_matrix_free(omp->reference_fyS[i]);
  }
  free(omp->reference_y);
  free(omp->reference_fy);
  free(omp->reference_yS);
  free(omp->reference_fyS);
 

  for (i=0;i<CT;i++){
    gsl_vector_free(omp->y[i]);
    gsl_vector_free(omp->fy[i]);
    gsl_matrix_free(omp->yS[i]);
    gsl_matrix_free(omp->fyS[i]);
    gsl_matrix_free(omp->oS[i]);

  }
  free(omp->y);
  free(omp->fy);
  free(omp->yS);
  free(omp->fyS);
  free(omp->oS);
  
  free(omp->data_row);
  free(omp->sd_data_row);
  free(omp->data);
  free(omp->sd_data);
  free(omp->input_u_row);
  free(omp->u);
  
  gsl_matrix_int_free(omp->norm_f);
  gsl_matrix_int_free(omp->norm_t);
  
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

  return EXIT_SUCCESS;
 
}
