#include "model_parameters_smmala.h"
/* Data will always be stored in blocks (matrices), there will be
 * vector views to identify specific rows.  Whether Data or data_block
 * is used in memory allocation is up to the data read functions.
 */
int init_E(ode_model_parameters *omp){
  int i;
  int C=omp->size->C;
  omp->E=(experiment**) malloc(sizeof(experiment*)*C);
  for (i=0;i<C;i++) {
    omp->E[i]=(experiment*) malloc(sizeof(experiment));
    omp->E[i]->t=NULL;
    omp->E[i]->init_y=NULL;
    omp->E[i]->input_u=NULL;
  }
  omp->ref_E=(experiment*) malloc(sizeof(experiment)); // just in case
  omp->ref_E->t=NULL;
  omp->ref_E->init_y=NULL;
  omp->ref_E->input_u=NULL;
  omp->Data=NULL;
  omp->sdData=NULL;
  for (i=0;i<C;i++) {
    omp->E[i]->data_block=NULL;
    omp->E[i]->sd_data_block=NULL;
    omp->E[i]->normalise=malloc(sizeof(normalisation_t));
    omp->E[i]->view=malloc(sizeof(view_t));    
  }
  omp->ref_E->data_block=NULL;
  omp->ref_E->sd_data_block=NULL;
  
  omp->S_approx=malloc(sizeof(sensitivity_approximation));
  
  omp->prior=malloc(sizeof(prior_t));
  printf("# %i+1 experiment structures allocated.\n",C);
  return GSL_SUCCESS;
}

int ode_model_parameters_link(ode_model_parameters *omp){
  int i,j,k=0;
  int T=0;
  int C=omp->size->C;
  int F=omp->size->F;

  for (i=0;i<C;i++){
    T=omp->E[i]->t->size; // there is one data block per experiment, with T rows
    //data
    if (omp->E[i]->data_block==NULL && omp->Data!=NULL){
      // submatrix views have to be created based on measurement times
      omp->E[i]->view->data_block=gsl_matrix_submatrix(omp->Data,k,0,T,F);
      omp->E[i]->data_block=&(omp->E[i]->view->data_block.matrix);
      // standard deviation; same treatment
      omp->E[i]->view->sd_data_block=gsl_matrix_submatrix(omp->sdData,k,0,T,F);
      omp->E[i]->sd_data_block=&(omp->E[i]->view->sd_data_block.matrix);
    }
    k+=T;
    // each row
    omp->E[i]->view->data_row=(gsl_vector_view*) malloc(sizeof(gsl_vector_view)*T);
    omp->E[i]->view->sd_data_row=(gsl_vector_view*) malloc(sizeof(gsl_vector_view)*T);
    omp->E[i]->data=(gsl_vector**) malloc(sizeof(gsl_vector*)*T);
    omp->E[i]->sd_data=(gsl_vector**) malloc(sizeof(gsl_vector*)*T);
    for (j=0;j<T;j++){
      omp->E[i]->view->data_row[j]=gsl_matrix_row(omp->E[i]->data_block,j);
      omp->E[i]->view->sd_data_row[j]=gsl_matrix_row(omp->E[i]->sd_data_block,j);
      omp->E[i]->data[j]=&(omp->E[i]->view->data_row[j].vector);
      omp->E[i]->sd_data[j]=&(omp->E[i]->view->sd_data_row[j].vector);
    }
  }
  //printf("data rows linked.\n");
  // shortcuts for input rows;
  return GSL_SUCCESS;
}

int ode_model_parameters_alloc(ode_model_parameters *omp){
  /* problem_size structure contains:
   * N: number of states
   * D: number of model parameters
   * F: number of output functions
   * U: number of input parameters
   * C: number of different experimental conditions
   * T: number of measured time points
   */
  int i,j;
  int N=omp->size->N; // number of state variables
  int C=omp->size->C; // number of experimental conditions
  int F=omp->size->F; // number of output functions
  int D=omp->size->D; // number of sampling parameters (to be estimated)
  int U=omp->size->U; // number of input parameters (known)
  int T=omp->size->T; // number of measurement time instances
  int P=omp->size->P; // number of total parameters D+U (a consistency check
	       // between ode_model and mcmc configuration file)
  
  printf("[model_parameters_alloc] allocating model parameters with: N=%i,\tP=%i,\tT=%i.\n",N,P,T);
  /* initialise gsl_vectors and matrices these will hold the
   * differential equation data for each experimental condition c and
   * time point t_j, like this: y[c*T+j]=gsl_vector(N)
   */
  omp->tmpF=gsl_vector_alloc(F);

  omp->S_approx->jacobian_y=gsl_matrix_alloc(N,N);
  omp->S_approx->jacobian_p=gsl_matrix_alloc(P,N);
  omp->S_approx->tau=gsl_vector_alloc(N);
  omp->S_approx->x=gsl_vector_alloc(N);
  omp->S_approx->r=gsl_vector_alloc(N);
  omp->S_approx->eJt=gsl_matrix_alloc(N,N);
  omp->S_approx->Jt=gsl_matrix_alloc(N,N);
  
  printf("[model_parameters_alloc] y, fy,yS, fyS, oS, yS0, nfy and nfyS.\n");
  for (i=0;i<C;i++){
    T=omp->E[i]->t->size;
    omp->E[i]->y=(gsl_vector **) malloc(sizeof(gsl_vector*)*T);
    omp->E[i]->fy=(gsl_vector **) malloc(sizeof(gsl_vector*)*T);
    omp->E[i]->yS=(gsl_matrix **) malloc(sizeof(gsl_matrix*)*T);
    omp->E[i]->fyS=(gsl_matrix **) malloc(sizeof(gsl_matrix*)*T);
    omp->E[i]->oS=(gsl_matrix **) malloc(sizeof(gsl_matrix*)*T);
    omp->E[i]->yS0=gsl_matrix_calloc(P,N);

    omp->E[i]->normalise->fy=malloc(sizeof(gsl_vector*)*T);
    omp->E[i]->normalise->fyS=malloc(sizeof(gsl_matrix*)*T);
    omp->E[i]->normalise->data=malloc(sizeof(gsl_vector*)*T);
    omp->E[i]->normalise->stdv=malloc(sizeof(gsl_vector*)*T);   
   
    for (j=0;j<T;j++){      
      omp->E[i]->y[j]=gsl_vector_alloc(N);
      omp->E[i]->fy[j]=gsl_vector_alloc(F);
      omp->E[i]->yS[j]=gsl_matrix_alloc(P,N);
      omp->E[i]->fyS[j]=gsl_matrix_alloc(P,F);
      omp->E[i]->oS[j]=gsl_matrix_alloc(D,F);
      // normalisation:
          omp->E[i]->normalise->fy[j]=gsl_vector_alloc(F);
	  omp->E[i]->normalise->fyS[j]=gsl_matrix_alloc(P,F);
	  omp->E[i]->normalise->data[j]=gsl_vector_alloc(F);
	  omp->E[i]->normalise->stdv[j]=gsl_vector_alloc(F);
    }
  }
  omp->ref_E->y=(gsl_vector **) malloc(sizeof(gsl_vector*)*T);
  omp->ref_E->fy=(gsl_vector **) malloc(sizeof(gsl_vector*)*T);
  omp->ref_E->yS=(gsl_matrix **) malloc(sizeof(gsl_vector*)*T);
  omp->ref_E->fyS=(gsl_matrix **) malloc(sizeof(gsl_vector*)*T);
  omp->ref_E->oS=(gsl_matrix **) malloc(sizeof(gsl_matrix*)*T);
  omp->ref_E->yS0=gsl_matrix_calloc(P,N);

  for (j=0;j<T;j++){
    omp->ref_E->y[j]=gsl_vector_alloc(N);
    omp->ref_E->fy[j]=gsl_vector_alloc(F);
    omp->ref_E->yS[j]=gsl_matrix_alloc(P,N);
    omp->ref_E->fyS[j]=gsl_matrix_alloc(P,F);
    omp->ref_E->oS[j]=gsl_matrix_alloc(D,F);
  }  
  /* during burn-in, we slowly increase beta from 0 to 1; if no
   * burn-in is performed, beta needs to be 1.0
   */
  omp->p=gsl_vector_alloc(P);

  
    // prior
  omp->prior->p=gsl_permutation_alloc((size_t) D);
  omp->prior->n=3;
  omp->prior->tmp=malloc(sizeof(gsl_vector*) * omp->prior->n);
  for (i=0;i< omp->prior->n;i++){
    omp->prior->tmp[i]=gsl_vector_alloc(D);
  }
  printf("[model_parameters_alloc] done.\n");
  return EXIT_SUCCESS;
}

int ode_model_parameters_free(ode_model_parameters *omp){
  int i,j,k;
  int C=omp->size->C;
  int T=omp->size->T;

  gsl_matrix_free(omp->S_approx->jacobian_y);
  gsl_matrix_free(omp->S_approx->jacobian_p);
  gsl_vector_free(omp->S_approx->tau);
  gsl_vector_free(omp->S_approx->x);
  gsl_vector_free(omp->S_approx->r);
  gsl_matrix_free(omp->S_approx->Jt);
  gsl_matrix_free(omp->S_approx->eJt);
  free(omp->S_approx);

  for (i=0;i<C;i++){
    free(omp->E[i]->view->data_row);
    free(omp->E[i]->view->sd_data_row);
    free(omp->E[i]->data);
    free(omp->E[i]->sd_data);
    for (j=0;j<T;j++){
      gsl_vector_free(omp->E[i]->y[j]);
      gsl_vector_free(omp->E[i]->fy[j]);
      gsl_matrix_free(omp->E[i]->yS[j]);
      gsl_matrix_free(omp->E[i]->fyS[j]);
      gsl_matrix_free(omp->E[i]->oS[j]);
    }
    free(omp->E[i]);
  }
  free(omp->E);
  for (j=0;j<T;j++){
    gsl_vector_free(omp->ref_E->y[j]);
    gsl_vector_free(omp->ref_E->fy[j]);
    gsl_matrix_free(omp->ref_E->yS[j]);
    gsl_matrix_free(omp->ref_E->fyS[j]);
    gsl_matrix_free(omp->ref_E->oS[j]);
  }
  free(omp->ref_E);
  gsl_vector_free(omp->p);
  free(omp->size);
  return EXIT_SUCCESS;
}
