#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_rng.h>
#include <string.h>
#include "../app/flatten.h"
#include <mpi.h>
#include <assert.h>
#include "smmala.h"
#include "../ode/ode_model.h"
#include "smmala_posterior.h"
#include "../app/diagnosis_output.h"
#include "model_parameters_smmala.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "omp.h"
int LogLikelihood(ode_model_parameters *mp, double *l, gsl_vector *dl, gsl_matrix *FI);

int ode_solver_process_sens(ode_solver *solver, double tout,
			    gsl_vector *y, gsl_vector* fy,
			    gsl_matrix *yS, gsl_matrix *fyS,
			    sensitivity_approximation *a){
  ode_model *model;
  model=solver->odeModel;
  
  if (ode_model_has_sens(model)){
    ode_solver_get_sens(solver, tout, yS->data);
  }  else {
    // printf("approximating sensitivity\n"); fflush(stdout);
    // approximate sensitivities using steady state assumptions;
    gsl_vector_view eJt_diag;
    //gsl_permutation *permutation;
    gsl_vector_view jacp_row;
    //int signum;
    int status;
    //int N=ode_model_getN(model);
    int P=ode_model_getP(model);
    int i;
    // get jacobian
    ode_solver_get_jac(solver,tout,y->data,fy->data,a->jacobian_y->data);
    ode_solver_get_jacp(solver,tout,y->data,fy->data,a->jacobian_p->data);
    gsl_matrix_transpose(a->jacobian_y); // now jacobian_y(i,j)=df[i]/df[j];
    //gsl_printf("Jac_y",a->jacobian_y,1);
    gsl_matrix_memcpy(a->Jt,a->jacobian_y);
    gsl_matrix_scale(a->Jt,tout);           // this is now Jacobian*t    
    // each row of jacp is a different parameter; row1 is df[i]/dp1
    status=gsl_linalg_QR_decomp(a->jacobian_y,a->tau);
    //    status=gsl_linalg_LU_decomp(jacobian_y, permutation, &signum);
    if (status!= GSL_SUCCESS){
      fprintf(stderr,"[sens approx] QR decomposition failed. %s\n",gsl_strerror(status));      
    } 
    for (i=0;i<P;i++){
      jacp_row=gsl_matrix_row(a->jacobian_p,i);
      // solve in place; jac_p will contain the solution
      status=gsl_linalg_QR_svx(a->jacobian_y, a->tau, &jacp_row.vector);
      if (status!=GSL_SUCCESS) {
	fprintf(stderr,"[sens approx] QR solution of linear equations failed: %s. Using minimum Norm solution.\n",gsl_strerror(status));
	if (gsl_linalg_QR_lssolve(a->jacobian_y, a->tau, &jacp_row.vector, a->x, a->r)!=GSL_SUCCESS) {
	  fprintf(stderr,"[sens approx] QR lssolve also failed. exiting.\n");
	  exit(-1);
	} else {
	  gsl_vector_memcpy(&jacp_row.vector,a->x);
	}	
      }
    }
    status=gsl_linalg_exponential_ss(a->Jt,a->eJt,GSL_PREC_SINGLE);
    if (status!=GSL_SUCCESS){
      // this is not yet considered stable by GSL :/
      fprintf(stderr,"[sens approx] matrix exponential failed. %s\n",gsl_strerror(status));
    }
    eJt_diag=gsl_matrix_diagonal(a->eJt);
    gsl_vector_add_constant(&eJt_diag.vector,-1.0);
    if (gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0,a->jacobian_p,a->eJt,1.0,yS)!=GSL_SUCCESS){
      fprintf(stderr,"[sens approx] dge Matrix·Matrix failed.\n");
    }
    //gsl_permutation_free(permutation);
  } // end if has sens
  
  if (ode_model_has_funcs_sens(model)){
    ode_solver_get_func_sens(solver, tout, y->data, yS->data, fyS->data);
  } // end if has func sens
  return GSL_SUCCESS;
}
//ode_solver_step(solver, gsl_vector_get(t,j), y, fy, yS, fyS, a);
int ode_solver_step(ode_solver *solver, double t, gsl_vector *y, gsl_vector* fy, gsl_matrix *yS, gsl_matrix *fyS, sensitivity_approximation *a){
  /* solves ODE for parameter vector p;
   * returns state vectors y
   */
  double tout;
  ode_model *model;
  model=solver->odeModel;
  int CVerror =  ode_solver_solve(solver, t, y->data, &tout);

  if (CVerror!=CV_SUCCESS) {
    fprintf(stderr, "ODE solver failed with ERROR = %i.\n",CVerror);
    // return a rejection message; the Likelihood is not defined for this argument;
    return GSL_EDOM;
  }
  // get sensitivities and output function values for the calculated ODE solutions
  if (ode_model_has_funcs(model)) {
    //printf("obtaining functions: "); fflush(stdout);
    ode_solver_get_func(solver, tout, y->data, fy->data);
    //gsl_vector_fprintf(stdout,fy,"%f, ");
    //printf("..done\n"); fflush(stdout);
  }
  else {
    fprintf(stderr,"ode model has no output functions");
    exit(-1);
  }
  //printf("obtaining sensitivities.\n"); fflush(stdout);
  ode_solver_process_sens(solver, tout, y, fy, yS, fyS, a);
  return GSL_SUCCESS;
}


/* Normalisation information: quantities are normalised in place:
 * NORMALISATION_TYPE: fy_i ← normalisation_function(fy_i,fy_j)
 * DATA_NORMALISED_BY_REFERENCE: fy_i(t_j,k,u)←fy_i(t_j,k,u)/fy_i(t_j,k,reference_u)
 * DATA_NORMALISED_BY_TIMEPOINT: fy_i(t_j,k,u)←fy_i(t_j,k,u)/fy_i(t_l,k,u) given l
 * DATA_NORMALISED_BY_STATE_VAR: fy_i(t_j,k,u)←fy_i(t_j,k,u)/fy_j(t_l,k,u),
 *                               given i, j and l
 * DATA_NORMALISED_INDIVIDUALLY: the model/data structure contains normalisation 
 *                               information in each experiment struct, so each 
 *                               experiment can be normalised in a different way. 
 *                               This is meant to supercede the other normalisation 
 *                               types.
 *
 *"reference" is not meant as in C/code jargon (pass by reference), it's more literal.
 * Note that all returned values/arrays from CVODES are with respect to
 * the model's parameters k=exp(x), while we require sensitivities
 * with respect to x, the mcmc variables. Therefore, chain rule factors appear in the 
 * below calculations: dk(x)/dx = k(x); df(k(x))/d(x) = df/dk dk/dx = df/dk * k
 */
/* this function finds the right normalisation vector and copies it
 * into the normalisation sub-structure of experiment E
 */
int get_normalising_vector_cpy(experiment *E, experiment *ref_E){
  int t,rt;
  int f,k;
  double fy_k;
  gsl_vector_view fyS_column, ref_column;
  int j=E->NormaliseByTimePoint;
  if(ref_E==NULL){
    // this function is called only if some normalisation is required,
    // so if ref_E is NULL; it defaults to E;
    ref_E=E;
  }
  int F=E->fy[0]->size;
  int T=E->t->size;
  int rT=ref_E->t->size;

  assert(j<T);

  if (j<0) assert(T==rT);
  // j<0 means that this experiment is not normalised by one specific
  // time point, but rather time-series-wise: fy[t]/rfy[t] for all t,
  // so fy amd rfy must be measured equally often T==rT. To make
  // sense, probably the measurements should occur at the same time;
  // but, we don't check that here. t is an index (SBtab timepoint).

  if (E->NormaliseByOutput==NULL){
    for (t=0;t<T;t++) {
      rt=j<0?t:j;
      gsl_vector_memcpy(E->normalise->fy[t],ref_E->fy[rt]);
      gsl_matrix_memcpy(E->normalise->fyS[t],ref_E->fyS[j]);
    }
  } else { // we have to re-order the outputs
    for (f=0;f<F;f++){
      k=gsl_vector_int_get(E->NormaliseByOutput,f);
      if (k<0 || k>=F) {
	for (t=0;t<T;t++) {
	  gsl_vector_set(E->normalise->fy[t],f,1.0);
	  fyS_column=gsl_matrix_column(E->normalise->fyS[t],f);
	  gsl_vector_set_zero(&(fyS_column.vector));
	}
      } else {
	for (t=0;t<T;t++) {
	  rt=j<0?t:j; // reference time
	  // copy the output values to normalisation struct
	  fy_k=gsl_vector_get(ref_E->fy[rt],k);
	  gsl_vector_set(E->normalise->fy[t],f,fy_k);
	  // copy sensitivity column to normalisation struct
	  fyS_column=gsl_matrix_column(E->normalise->fyS[t],f);
	  ref_column=gsl_matrix_column(ref_E->fyS[rt],k);
	  gsl_vector_memcpy(&(fyS_column.vector),&(ref_column.vector));
	}
      }      
    }
  }  
  return GSL_SUCCESS;
}

/* A more general normalisation routine. It reads out normalisation
 * information from the actual experiment struct and treats every
 * experiment differently, as specified. (work in progress)
 */
int normalise(ode_model_parameters *mp){
  int c,j,k;
  int C=mp->size->C;
  int T=mp->size->T;
  int D=mp->size->D;
  int F=mp->size->F;
  assert(T>0);
  gsl_vector *fy, *rfy;
  gsl_matrix *fyS, *oS, *rfyS;
  gsl_matrix_view fyS_D_sub, rfyS_D_sub; // the first D rows of fyS;
  gsl_matrix *fySD, *rfySD; // pointer to the above;
  gsl_vector *v;
  gsl_matrix *M;
  gsl_vector_view M_row;
  /* gsl_vector_view fyS_k_row; //  */
  /* gsl_vector *fyS_k; // the raw fy sensitivity for one parameter k */
  /* gsl_vector_view rfyS_k_row; //  */
  /* gsl_vector *rfyS_k; // the raw fy sensitivity for one parameter k */
  gsl_vector_view oS_k_row; // 
  gsl_vector *oS_k; // the normalised output sensitivity for one parameter k
  int i,nt;
  //int rt;
  experiment *ref_E=NULL;
  v=mp->tmpF;
  M=mp->tmpDF;
  assert(v && M);
  //printf("[normalise] normalising all Experiments according to the flags: NormaliseByExperiment, NormaliseByTimePoint and NormaliseByOutput.\n");
  
  for (c=0;c<C;c++){
    nt=mp->E[c]->t->size;
    assert(nt<=T);
    if (NEEDS_NORMALISATION(mp->E[c])){
      i=mp->E[c]->NormaliseByExperiment;
      ref_E=(i<0 || i>C)?NULL:(mp->E[i]);
      //rt=mp->E[c]->NormaliseByTimePoint;
      // the result of get_normalising_vector_cpy wil be stored in the
      // experiment's normalisation struct
      get_normalising_vector_cpy(mp->E[c], ref_E);
      for (j=0;j<nt;j++){
	fy=mp->E[c]->fy[j];
	fyS=mp->E[c]->fyS[j];
	oS=mp->E[c]->oS[j];
	rfy=mp->E[c]->normalise->fy[j];
	if (!gsl_vector_ispos(rfy)){
	  gsl_vector_add_constant(rfy,1E-10);
	}
	rfyS=mp->E[c]->normalise->fyS[j];
	gsl_vector_div(fy,rfy);
	// here, we create matrix views of the right size for
	// convenience: fyS is of size P×F, but for mcmc purposes we
	// don't need input sensitivities, and p=[x,u]: P=D+U (all
	// parameters = "mcmc related" parameters x and input parameters u)
	fyS_D_sub=gsl_matrix_submatrix(fyS,0,0,D,F);
	rfyS_D_sub=gsl_matrix_submatrix(rfyS,0,0,D,F);
	fySD=&(fyS_D_sub.matrix);
	rfySD=&(rfyS_D_sub.matrix);
	gsl_matrix_memcpy(oS,fySD);
	gsl_matrix_memcpy(M,rfySD);
	for (k=0;k<D;k++) {
	  M_row=gsl_matrix_row(M,k);
	  gsl_vector_mul(&(M_row.vector),fy);
	}
	gsl_matrix_sub(oS,M);
	for (k=0;k<D;k++) {
	  oS_k_row=gsl_matrix_row(oS,k);
	  oS_k=&(oS_k_row.vector);
	  gsl_vector_div(oS_k,rfy);
	  gsl_vector_scale(oS_k,gsl_vector_get(mp->p,k)); // because we are in log-space
	}
      }	     
    } else { // so, normalisation is not needed
      for (j=0;j<nt;j++){
	fy=mp->E[c]->fy[j];
	fyS=mp->E[c]->fyS[j];
	oS=mp->E[c]->oS[j];
	fyS_D_sub=gsl_matrix_submatrix(fyS,0,0,D,F);
	gsl_matrix_memcpy(oS,&(fyS_D_sub.matrix));
	for (k=0;k<D;k++){
	  oS_k_row=gsl_matrix_row(oS,k);
	  oS_k=&(oS_k_row.vector);	  
	  gsl_vector_scale(oS_k,gsl_vector_get(mp->p,k)); // because we are in log-space
	}
      }
    }
  }
  return GSL_SUCCESS;
}

int assign_oS_fyS(ode_model_parameters *mp){
  int j,k,c;
  int C,T;
  int D,F;
  D=mp->size->D;
  F=mp->size->F;
  C=mp->size->C;
  gsl_vector *t;
  gsl_matrix_view fyS_sub;
  gsl_vector_view oS_row;
  //printf("assigning oS←fyS\n");
  for (c=0;c<C;c++){
    t=mp->E[c]->t;
    T=t->size;
    for (j=0;j<T;j++){
      fyS_sub=gsl_matrix_submatrix(mp->E[c]->fyS[j],0,0,D,F);
      gsl_matrix_memcpy(mp->E[c]->oS[j],&(fyS_sub.matrix));
      for (k=0;k<D;k++){
	oS_row=gsl_matrix_row(mp->E[c]->oS[j],k);
	gsl_vector_scale(&oS_row.vector,gsl_vector_get(mp->p,k));
      }
    }
  }
  return GSL_SUCCESS;
}


int save_simulation_results(ode_model_parameters *omp){
  int c;
  int C=omp->size->C;
  herr_t EC=0;
  hid_t file_id=H5Fcreate("SimulationResult.h5",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  ode_model *model=omp->solver[0]->odeModel;
  
  assert(file_id>0);
  printf("[save_simulation_results] hdf5 file opened\n");
  const char **y_name=ode_model_get_var_names(model);
  const char **p_name=ode_model_get_param_names(model);
  const char **f_name=ode_model_get_func_names(model);

  char *y_names=flatten(y_name, (size_t) omp->size->N, "; ");
  char *p_names=flatten(p_name, (size_t) omp->size->P, "; ");
  char *f_names=flatten(f_name, (size_t) omp->size->F, "; ");

  herr_t NameWriteError=0;
  NameWriteError&=H5LTmake_dataset_string(file_id,"StateVariableNames",y_names);
  NameWriteError&=H5LTmake_dataset_string(file_id,"ParameterNames",p_names);
  NameWriteError&=H5LTmake_dataset_string(file_id,"OutputFunctionNames",f_names);
  assert(NameWriteError==0);
  free(y_names);
  free(p_names);
  free(f_names);
  
  hsize_t t_size[1]={1}; // one measurement time index
  hsize_t y_size[1]={omp->size->N};
  hsize_t u_size[1]={omp->size->U};
  hsize_t p_size[1]={omp->size->D};
  hsize_t fy_size[1]={omp->size->F};
  hsize_t yS_size[2]={omp->size->P,omp->size->N};
  hsize_t fyS_size[2]={omp->size->P,omp->size->F};
  hsize_t oS_size[2]={omp->size->D,omp->size->F};
  hsize_t data_size[2];

  int t,T=0;
  int F=omp->size->F;

  hid_t experiment_id, time_id;
  char ExperimentName[128];
  char TimeName[128];
  size_t n=127,nout=-1;
  printf("[save_simulation_results] There are %i simulations to save.\n",C);
  EC&=H5LTmake_dataset_double(file_id,"ModelParameters", 1, p_size, omp->p->data);
  assert(EC==0);

  for (c=0;c<C;c++){
    model=omp->solver[c]->odeModel;
    T=omp->E[c]->t->size;    
    t_size[0]=T;
    data_size[0]=T;
    data_size[1]=F;
    nout=snprintf(ExperimentName,n,"/Experiment%i",c);
    assert(nout>0 && nout<n);
    experiment_id=H5Gcreate2(file_id,ExperimentName, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    EC&=H5LTmake_dataset_double(experiment_id,"time", 1, t_size, omp->E[c]->t->data);
    assert(EC==0);

    EC&=H5LTmake_dataset_double(experiment_id,"data_block", 2, data_size, omp->E[c]->data_block->data);
    assert(EC==0);
    EC&=H5LTmake_dataset_double(experiment_id,"stdv_block", 2, data_size, omp->E[c]->sd_data_block->data);
    assert(EC==0);      
    EC&=H5LTmake_dataset_double(experiment_id,"InputParameters", 1, u_size, omp->E[c]->input_u->data);
    assert(EC==0);
    t_size[0]=1;
    data_size[0]=1;
    data_size[1]=F;
    for (t=0;t<T;t++){
      nout=snprintf(TimeName,n,"TimePoint%i",t);
      assert(nout>0 && nout<n);
      time_id=H5Gcreate2(experiment_id,TimeName, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

      EC&=H5LTmake_dataset_double(time_id,"data_row", 2, data_size, omp->E[c]->data[t]->data);
      assert(EC==0);
      EC&=H5LTmake_dataset_double(time_id,"stdv_row", 2, data_size, omp->E[c]->sd_data[t]->data);
      assert(EC==0);      
      EC&=H5LTmake_dataset_double(time_id,"MeasurementTime", 1, t_size, &(omp->E[c]->t->data[t]));
      assert(EC==0);
      EC&=H5LTmake_dataset_double(time_id,"StateVariables", 1, y_size, omp->E[c]->y[t]->data);
      assert(EC==0);
      EC&=H5LTmake_dataset_double(time_id,"OutputFunctions", 1, fy_size, omp->E[c]->fy[t]->data);
      assert(EC==0);
      EC&=H5LTmake_dataset_double(time_id,"StateSensitivityRaw", 2, yS_size, omp->E[c]->yS[t]->data);
      assert(EC==0);
      EC&=H5LTmake_dataset_double(time_id,"FunctionSensitivityRaw", 2, fyS_size, omp->E[c]->fyS[t]->data);
      assert(EC==0);
      EC&=H5LTmake_dataset_double(time_id,"FunctionSensitivityLogSpace", 2, oS_size, omp->E[c]->oS[t]->data);
      assert(EC==0);      
      EC&=H5Gclose(time_id);
    }
    EC&=H5Gclose(experiment_id);
  }
  EC&=H5Fclose(file_id);
  assert(EC==0);
  return EC;
}


int LogLikelihood(ode_model_parameters *mp, double *l, gsl_vector *grad_l, gsl_matrix *fisher_information){

  /* Here, the ode integration is done
   * 
   * «l» is a scalar, the return slot of the log-likelihood
   * dl is the parameter gradient of the log likelihood (size=D)
   *
   * cvode returns functions «fy» of the state variables «y».
   * 
   * Note P is the cvode related number of ode parameters
   * while D is the number of MCMC related parameters
   * P=D+U
   */
  int i,j,c,T;
  size_t k,K;
  double l_t_j;
  double e_t; //event_time
  int status=GSL_SUCCESS;
  gsl_vector *y,*fy;
  gsl_matrix *yS,*fyS;
  gsl_vector *t;
  gsl_vector_view oS_row;
  gsl_vector_view input_part;
  int i_flag=GSL_SUCCESS; // ODE (i)ntegration error flag
  ode_model *model;
  ode_solver **solver;
  sensitivity_approximation *a;
  solver=mp->solver;
  assert(solver);
  assert(mp->size);
  int N=get_number_of_state_variables(mp);
  int D=get_number_of_MCMC_variables(mp);
  int P=get_number_of_model_parameters(mp);
  //int F=get_number_of_model_outputs(mp);
  int U=get_number_of_model_inputs(mp);
  int C=get_number_of_experimental_conditions(mp);
  assert(P>0 && N>0);
  /* The ODE model integration below takes by far the most time to
   * calculate.  This is the only bit of the code that needs to be
   * parallel apart from parallel tempering done by mpi.
   */
  //printf("[%s] P=%i, U=%i.\n",__func__,P,U);
#pragma omp parallel for private(model,j,e_t,k,K,a,y,fy,yS,fyS,t,T,input_part) reduction(|:i_flag)
  for (c=0; c<C; c++){/* loop over different experimental conditions */
    model=solver[c]->odeModel;    
    a=mp->S_approx[c];
    assert(mp->E[c]->p);
    gsl_vector_memcpy(mp->E[c]->p,mp->p); // events may modify this vector
    /* each loop iteration gets a different input vector, because an
       event may change the input parameters.
     */
    input_part=gsl_vector_subvector(mp->E[c]->p,D,U);    
    gsl_vector_memcpy(&(input_part.vector),mp->E[c]->input_u);
    ode_solver_reinit(solver[c], mp->t0, mp->E[c]->init_y->data, N,
		      mp->E[c]->p->data,
		      mp->E[c]->p->size);
    if (ode_model_has_sens(model)){
      ode_solver_reinit_sens(solver[c], mp->E[c]->yS0->data, P, N);
    }
    t=mp->E[c]->t;
    T=t->size;
    for (j=0; j<T; j++){
      y=mp->E[c]->y[j]; //printf("y: %i, %zi\n",mp->size->N,y->size);
      fy=mp->E[c]->fy[j]; //printf("fy: %i, %zi\n",mp->size->F,fy->size);
      yS=mp->E[c]->yS[j]; //printf("yS: %i, %zi\n",mp->size->N*P,yS->size1*yS->size2);
      fyS=mp->E[c]->fyS[j]; //printf("fyS: %i, %zi\n",mp->size->F*P,fyS->size1*fyS->size2);
      /* 1. process all events that precede t(j) */
      assert(y && fy && yS && fyS);
      if (mp->E[c]->event_list && mp->E[c]->before_t && mp->E[c]->before_t[j]){
	K=mp->E[c]->before_t[j]->size; // number of events prior to t[j]
	for (k=0;k<K;k++){
	  e_t=mp->E[c]->before_t[j]->event[k]->t;
	  i_flag|=ode_solver_step(solver[c], e_t, y, fy, yS, fyS, a);
	  event_apply(mp->E[c]->before_t[j]->event[k],y,mp->E[c]->p,yS);
	}
      }
      /* 2. advance the state to the measurement time t(j) */
      i_flag|=ode_solver_step(solver[c], gsl_vector_get(t,j), y, fy, yS, fyS, a);
    }
  }
  
  if (i_flag!=GSL_SUCCESS){
    l[0]=-INFINITY;
    return i_flag;
  }else{
    //printf("pre normalisation\n");
    //printf_omp(mp);    
    switch (mp->normalisation_type){
    case DATA_NORMALISED_INDIVIDUALLY:
      normalise(mp);
      break;      
    case DATA_IS_ABSOLUTE:
      assign_oS_fyS(mp);
      break;
    default:
      fprintf(stderr,"unknown normalisation method: %i\n",mp->normalisation_type);
      exit(-1);
    }
    //printf("post normalisation\n");
    //printf_omp(mp);
    // save_simulation_results(mp); fflush(stdout);
    //initialise all return values
    // log-likelihood
    l[0]=mp->pdf_lognorm;    
    // gradient of the log-likelihood
    gsl_vector_set_zero(grad_l);
    gsl_matrix_set_zero(fisher_information);  
  
    for (c=0; c<C; c++){
      t=mp->E[c]->t;
      T=t->size;
      if (mp->E[c]->lflag){
	for (j=0; j<T; j++){
	  /* Calculate log-likelihood value
	   */      
	  gsl_vector_sub(mp->E[c]->fy[j],mp->E[c]->data[j]);
	  gsl_vector_div(mp->E[c]->fy[j],mp->E[c]->sd_data[j]); // (fy-data)/sd_data
	  l_t_j=0; // l_t_j =?= mp->E[c]->fy[j]->size * (0.5*log(2*M_PI) + gsl_vector_prod(mp->E[c]->sd_data[j]));
	  if (gsl_blas_ddot(mp->E[c]->fy[j],mp->E[c]->fy[j], &l_t_j)!=GSL_SUCCESS){
	    printf("[%s] blas ddot was unsuccessful\n",__func__);
	    exit(-1);
	  } // sum((fy-data)²/sd_data²)
	  if (isnan(l_t_j)){
	    fprintf(stderr,"[%s] likelihood evaluation resulted in nan.\n",__func__);
	    printf_omp(mp);
	    fflush(stdout);
	    fflush(stderr);
	    exit(-1);
	  }
	  l[0]+=-0.5*l_t_j;
	  /* Calculate The Likelihood Gradient and Fisher Information:
	   */
	  gsl_vector_div(mp->E[c]->fy[j],mp->E[c]->sd_data[j]); // (fy-data)/sd_data²
	  //gsl_dgemv(TransA,alpha,A,x,beta,y)
	  status=gsl_blas_dgemv(CblasNoTrans,-1.0,mp->E[c]->oS[j],mp->E[c]->fy[j],1.0,grad_l);
	  // grad_l += - (fy-data)/sd_data² * dfy/dtheta
	  if (status!=GSL_SUCCESS){
	    gsl_printf("oS",mp->E[c]->oS[j],1);
	    gsl_printf("(fy-data)/sd_data²",mp->E[c]->fy[j],0);
	    fprintf(stderr,"[%s] dgemv was unsuccessful: %i %s\n",__func__,status,gsl_strerror(status));
	    //exit(-1);
	  }
	  /* Calculate the Fisher information
	   */
	  for (i=0;i<D;i++) {
	    oS_row=gsl_matrix_row(mp->E[c]->oS[j],i);
	    assert(gsl_vector_ispos(mp->E[c]->sd_data[j]));
	    gsl_vector_div(&(oS_row.vector),mp->E[c]->sd_data[j]);
	  }
	  gsl_blas_dgemm(CblasNoTrans,
			 CblasTrans, 1.0,
			 mp->E[c]->oS[j],
			 mp->E[c]->oS[j], 1.0,
			 fisher_information);	
	} // end for loop for time points
      } // if lflag
    } //end for different experimental conditions (i.e. inputs)
  }
  //printf("[LogLikelihood] done.\n"); fflush(stdout);
  return i_flag & status;
}

int display_prior_information(const void *Prior_str){
  prior_t *prior=(prior_t*) Prior_str;
  int D=prior->mu->size;
  //gsl_vector *prior_diff;
  assert(prior!=NULL);
  // get the ode model's prior type
  int opt=prior->type;
  printf("[prior info] prior structure of type %i:",opt);  
  if (PTYPE(opt,PRIOR_IS_GAUSSIAN)){
    printf("Gaussian ");
    if (PTYPE(opt,PRIOR_IS_MULTIVARIATE)){
      printf("multivariate ");
      // 1. dprior=-Sigma\(x-mu)
      if (PTYPE(opt,PRIOR_SIGMA_GIVEN)){
	printf("with given Sigma, LU decomposed: \n");
	// solve linear equation (Ax=b <-> x=A\b):
	assert(prior->Sigma_LU!=NULL && prior->p!=NULL);
	gsl_printf("Sigma LU",prior->Sigma_LU,GSL_IS_DOUBLE | GSL_IS_MATRIX);
	gsl_permutation_fprintf(stdout, prior->p, " %i ");
      } else if (PTYPE(opt,PRIOR_PRECISION_GIVEN)) {
	// do (matrix * vector) operation
	printf("with given inverted Sigma: \n");
	assert(prior->inv_cov!=NULL && prior->inv_cov->size1==D);
	gsl_printf("inv(Sigma)",prior->inv_cov,GSL_IS_DOUBLE | GSL_IS_MATRIX);
      } else {
	fprintf(stderr,"[LogPrior] type %i is GAUSSIAN but neither SIGMA nor PRECISION is given.\n",opt);
      }
    } else {
      // alternatively: product distribution ...
      printf("product distribution:\n");
      gsl_printf("sigma",prior->sigma,GSL_IS_DOUBLE | GSL_IS_VECTOR);
    }
  } else{
    printf("[prior info] prior is unknown.\n");
    fprintf(stderr,"[prior] type %i is not handled right now.\n",prior->type);
  }
  gsl_printf("mu",prior->mu,GSL_IS_DOUBLE | GSL_IS_VECTOR);
  fflush(stdout);
  return GSL_SUCCESS;
}

int LogPrior(const prior_t *prior, const gsl_vector *x, double *prior_value, gsl_vector *dprior, gsl_matrix *fi){
  assert(prior_value!=NULL && dprior!=NULL);
  int D=prior->mu->size;
  gsl_vector *prior_diff;
  assert(prior!=NULL);
  prior_value[0]=0;                     // the prior distribution related sum of squares
  prior_diff=prior->tmp[0];
  assert(prior_diff!=NULL && prior->n>0);
  // 0. prior_diff = (x-mu)
  gsl_vector_memcpy(prior_diff,x);
  gsl_vector_sub(prior_diff,prior->mu);
  //printf("[LogPrior] calculating prior.\n");
  int gsl_status=GSL_SUCCESS;
  // get the ode model's prior type
  int opt=prior->type;
  gsl_matrix_set_zero(fi);
  gsl_vector_set_zero(dprior); // prior gradient
  gsl_vector_view fi_diag;

  // useful function to incorporate later on:
  // include "gsl/gsl_randist.h"
  // double gsl_ran_gaussian_pdf(double x, double sigma);
  //
  
  if (PTYPE(opt,PRIOR_IS_GAUSSIAN)){
    if (PTYPE(opt,PRIOR_IS_MULTIVARIATE)){
      // 1. dprior=-Sigma\(x-mu)
      if (PTYPE(opt,PRIOR_SIGMA_GIVEN)){
	// solve linear equation (Ax=b <-> x=A\b):
	assert(prior->Sigma_LU!=NULL && prior->p!=NULL);
	gsl_status &= gsl_linalg_LU_solve(prior->Sigma_LU, prior->p, prior_diff,dprior);
	gsl_vector_scale(dprior,-1.0);
	// fisher information.... Sigma\eye(D);
	gsl_status &= gsl_linalg_LU_invert(prior->Sigma_LU, prior->p,fi);
      } else if (PTYPE(opt,PRIOR_PRECISION_GIVEN)) {
	// do (matrix * vector) operation
	assert(prior->inv_cov!=NULL && prior->inv_cov->size1==D);
	gsl_status &= gsl_blas_dsymv(CblasUpper,-1.0,prior->inv_cov,prior_diff,0.0,dprior);
	gsl_status &= gsl_matrix_memcpy(fi,prior->inv_cov);
      } //else {
	//fprintf(stderr,"[prior] type %i is GAUSSIAN but neither SIGMA nor PRECISION is given.\n",opt);
      //}
    } else {
      // alternatively: product distribution ... Gaussians
      gsl_status &= gsl_vector_memcpy(dprior,prior_diff);     // (x-mu)
      gsl_status &= gsl_vector_scale(dprior,-1.0);            // -(x-mu)
      gsl_status &= gsl_vector_div(dprior,prior->sigma); // -(x-mu)/sigma
      gsl_status &= gsl_vector_div(dprior,prior->sigma); // -(x-mu)/sigma²
      // fisher information
      fi_diag=gsl_matrix_diagonal(fi);
              gsl_matrix_set_identity(fi);
      gsl_status &= gsl_vector_div(&(fi_diag.vector),prior->sigma);
      gsl_status &= gsl_vector_div(&(fi_diag.vector),prior->sigma);
    }
    // 2. log_prior_value = 0.5 * [(-1)(x-mu)*Sigma\(x-mu)]
    gsl_status &= gsl_blas_ddot(prior_diff,dprior,prior_value);
    prior_value[0]*=0.5;
    assert(gsl_status==GSL_SUCCESS);
  } else{
    fprintf(stderr,"[prior] type %i is not handled right now.\n",prior->type);
    exit(-1);
  }
  //printf("[LogPrior] (%g) done.\n",prior_value[0]);
  //fflush(stdout);
  return gsl_status;
}

/* Calculates the unormalised log-posterior fx, the gradient dfx, the
 * Fisher information FI for Markov chain position x.
 */
int LogPosterior(const double beta, const gsl_vector *x,  void* model_params, double *fx, gsl_vector **dfx, gsl_matrix **FI){
  ode_model_parameters* omp =  (ode_model_parameters*) model_params;
  int i;
  int status=GSL_SUCCESS;
  int D=omp->prior->mu->size;
  double x_i;
  double p_i;
  assert(D>0);
  //gsl_sf_result result;
  //int gsl_status;
  for (i=0;i<D;i++) {
    x_i=gsl_vector_get(x,i);
    //gsl_status=gsl_sf_exp_e(x_i, &result);
    p_i=gsl_sf_exp(x_i);
    /* if (gsl_status==GSL_SUCCESS && result.val>result.err){ */
    /*   p_i=result.val; */
    /* }else{ */
    /*   fprintf(stderr,"[LogPosterior] gsl_sf_exp(x[%i])=%f±%f; failed: x=%g.\n",i,result.val,result.err,x_i); */
    /*   p_i=0.0; */
    /* } */
    gsl_vector_set(omp->p,i,p_i);
  }
  status &= LogLikelihood(omp, &fx[i_likelihood], dfx[i_likelihood], FI[i_likelihood]); 
  status &= LogPrior(omp->prior, x, &fx[i_prior], dfx[i_prior], FI[i_prior]);
  //printf("PosteriorFI: prior_value=%f\n",prior_value);

  fx[i_posterior]=fx[i_likelihood];
  fx[i_posterior]*=beta;
  fx[i_posterior]+=fx[i_prior];

  gsl_vector_memcpy(dfx[i_posterior],dfx[i_likelihood]);
  gsl_vector_scale(dfx[i_posterior],beta);
  gsl_vector_add(dfx[i_posterior],dfx[i_prior]);

  gsl_matrix_memcpy(FI[i_posterior],FI[i_likelihood]);
  gsl_matrix_scale(FI[i_posterior],beta*beta);
  gsl_matrix_add(FI[i_posterior],FI[i_prior]);
  //printf("[LogPosterior] done.\n"); fflush(stdout);
  return status;
}


