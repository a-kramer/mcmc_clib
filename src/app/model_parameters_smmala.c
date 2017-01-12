#include "model_parameters_smmala.h"

int ode_model_experiment_link(experiment *E){
  int i;
  int T=E->t->size;
  E->data_row=(gsl_vector_view*) malloc(sizeof(gsl_vector_view)*T);
  E->sd_data_row=(gsl_vector_view*) malloc(sizeof(gsl_vector_view)*T);
  // list of vector shorcuts to the views, no need to allocate the
  // vectors themselves
  omp->data=(gsl_vector**) malloc(sizeof(gsl_vector*)*T);
  omp->sd_data=(gsl_vector**) malloc(sizeof(gsl_vector*)*T);
  //printf("assign Data row views\n");

  for (i=0;i<T;i++){
    omp->data_row[i]=gsl_matrix_row(omp->Data,i);
    omp->data[i]=&(omp->data_row[i].vector);
    omp->sd_data_row[i]=gsl_matrix_row(omp->sdData,i);
    omp->sd_data[i]=&(omp->sd_data_row[i].vector);
  }
  // shortcuts for input rows;
  return GSL_SUCCESS;

}

int ode_model_parameters_alloc(ode_model_parameters *omp, const problem_size *ps){
  /* problem_size structure contains:
   * N: number of states
   * D: number of model parameters
   * F: number of output functions
   * U: number of input parameters
   * C: number of different experimental conditions
   * T: number of measured time points
   */
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
  omp->E=(experiment*) malloc(sizeof(experiment)*C);

  /* during burn-in, we slowly increase beta from 0 to 1; if no
   * burn-in is performed, beta needs to be 1.0
   */
  omp->beta=1.0;
  
  omp->p=gsl_vector_alloc(P);
  // time
  printf("gsl alloc t\n");
  omp->t=gsl_vector_alloc(T);

  // Data (t is the fast index, c the slow one). One huge block of Data:
  //printf("gsl alloc Data\n");

  // prior 
  omp->prior_tmp_a=gsl_vector_alloc(D);
  omp->prior_tmp_b=gsl_vector_alloc(D);
  omp->prior_mu=gsl_vector_alloc(D);
  omp->prior_inverse_cov=gsl_matrix_alloc(D,D);

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

  gsl_matrix_free(omp->initial_conditions_y);
  free(omp->init_y_view);
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
