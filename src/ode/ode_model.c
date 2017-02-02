/*
 *  ode_model.c
 *  odeSolve
 *
 *  Created by Vassilios Stathopoulos on 28/10/2011.
 *  Copyright 2011 Computing Science. All rights reserved.
 *
 */

#include <dlfcn.h>
#include <libgen.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "ode_model.h"


/* Auxiliary function for getting the symbol name. 
 * Implementation of this function is subject to change.
 */
static char* get_vf_name(const char* filename){
  char* _post = "_odeModel";
  char* lib_name = basename(filename);
  size_t len = strlen(lib_name);
  size_t i;
  /* stop once a '_c' or '.' is found. Take that as the model name */
  for (i = 0; i < len; i++)
    if( ((lib_name[i] == '_') && (lib_name[i+1] == 'c')) || (lib_name[i] == '.') )
      break;
  
  size_t new_size = i + strlen(_post) + 1;
  char* ret = (char*) malloc( new_size * sizeof(char) );
  strncat(ret, lib_name, i);
  strcat(ret, _post);
  return ret;
}

ode_model* ode_model_loadFromFile(const char *filename){
	
  void* odeLibrary = dlopen(filename, RTLD_LOCAL | RTLD_LAZY);				/* resource acc */
  
  if (odeLibrary == 0){
    fprintf(stderr, "Library %s could not be loaded.\n",filename);
    char * err = dlerror();
    if(err!=0)
      fprintf(stderr, "%s",err);
    /* exit(1); */
    return 0;
  }
  
  /* get lib name */
  char* model_name = get_vf_name(filename);
  ode_model* odeModel = (ode_model*)dlsym(odeLibrary, model_name);
  free(model_name);
  
  if(odeModel==0){
    fprintf(stderr, "Could not find ode_model symbol in library %s.\n", filename);
    char * err = dlerror();
    if(err!=0)
      fprintf(stderr, "%s",err);
    
    dlclose(odeLibrary);													/* resource free */
    /* exit(1); */
    return 0;
  }
  
  /* store a pointer to the loaded library so we can close it later with dlclose */
  odeModel->dylib = odeLibrary;
  return odeModel;
	
}

double ode_model_getN(const ode_model* model){
  return model->N;
}

double ode_model_getP(const ode_model* model){
  return model->P;
}

double ode_model_getF(const ode_model* model){
  return model->F;
}

int ode_model_has_sens(const ode_model* model){
  return (model->vf_sens)? 1 : 0 ;
}	

int ode_model_has_funcs(const ode_model* model){
  return (model->vf_func)? 1 : 0 ;
}

int ode_model_has_funcs_sens(const ode_model* model){
  return (model->vf_func_sens)? 1 : 0 ;
}

const char** ode_model_get_var_names(const ode_model* model){
  return (const char**) model->v_names;
}

const char** ode_model_get_param_names(const ode_model* model){
  return (const char**) model->p_names;
}

const char** ode_model_get_func_names(const ode_model* model){
  return (const char**) model->f_names;
}



void ode_model_get_initial_conditions(const ode_model* model, double* y0, int lenY0){
  int i,N;
  N = model->N;
  N = (lenY0 < N)? lenY0 : N ;
  
  double* m_v = model->v;
  for (i=0; i<N; i++) {
    y0[i] = m_v[i];
  }
}

void ode_model_get_default_params(const ode_model* model, double* params, int lenParams){
  int i,P;
  P = model->P;
  P = (lenParams < P)? lenParams : P ;
  
  double* m_p = model->p;
  for (i=0; i<P; i++) {
    params[i] = m_p[i];
  }
}



void ode_model_free(ode_model* model){
  if(model->dylib != 0){
    dlclose(model->dylib);
  }
}


ode_solver*	ode_solver_alloc(ode_model* model){
  
  ode_solver* solver = (ode_solver*) malloc( sizeof(ode_solver) );				/* alloc */
  if (solver == 0){
    /* TODO: write a proper error handler */
    fprintf(stderr,"malloc failed to allocate memory for ode_solver\n");
    return 0;
  }
  
  solver->cvode_mem = CVodeCreate(CV_BDF,CV_NEWTON);								/* alloc */
  if( solver->cvode_mem == 0){
    /* TODO: write a proper error handler */
    fprintf(stderr,"CVodeCreate failed to allocate memory in ode_solver for cvode_mem.\n");
    
    free(solver);
    return 0;
  }
  
  int P = ode_model_getP(model);
  solver->params = (double*) malloc( sizeof(double) * P );						/* alloc */
  if( solver->params == 0 ){
    /* TODO: write a proper error handler */
    fprintf(stderr,"malloc failed to allocate memory in ode_solver for params.\n");
    CVodeFree(solver->cvode_mem);
    free(solver);
    return 0;
  }
  ode_model_get_default_params(model, solver->params, P);	
  
  int N = ode_model_getN(model);
  solver->odeModel = model;
  solver->y = N_VNewEmpty_Serial(N);												/* alloc */
  NV_DATA_S(solver->y) = solver->odeModel->v;
  solver->yS = 0;
  
  return solver;
}

void ode_solver_setErrTol(ode_solver* solver, const double rel_tol,  double* abs_tol, const int abs_tol_len){
  
  if( (abs_tol_len != 1) && (abs_tol_len != solver->odeModel->N)){
    fprintf(stderr,"ode_solver_setErrTol: length of abs_tol must be 1 or equal to the number of variables in the ode model.\n");
    return ;
  }
  
  /* set tollerances to the cvode_mem internal structure */
  if ( abs_tol_len == 1 )
    CVodeSStolerances(solver->cvode_mem, rel_tol, abs_tol[0]);
  else{
    N_Vector abs_tol_vec =  N_VNewEmpty_Serial(abs_tol_len);		/* alloc */
    NV_DATA_S(abs_tol_vec) = abs_tol;
    CVodeSVtolerances(solver->cvode_mem, rel_tol, abs_tol_vec);
    N_VDestroy_Serial(abs_tol_vec);									/* free */
  }
}


void ode_solver_init(ode_solver* solver, const double t0,  double* y0, int lenY, double* p, int lenP ){
  int i,flag;
  /* Get parameters */
  
  if (p != 0){
    if (lenP != solver->odeModel->P) {
      fprintf(stderr,"ode_solver_init: lenP must be equal %d, the number of parameters in the ode model.\n",solver->odeModel->P);
      return ;
    }
    
    for(i = 0; i < solver->odeModel->P; i++)
      solver->params[i] = p[i];
  }
  
  /* Get initial conditions */
  if(y0 != 0){
    if( lenY != solver->odeModel->N ){
      fprintf(stderr,"ode_solver_init: lenY must be equal %d, the number of variables in the ode model.\n",solver->odeModel->N);
      return ;
    }
    NV_DATA_S(solver->y) = y0;
  }
  
  /* initialise */
  flag = CVodeInit(solver->cvode_mem, solver->odeModel->vf_eval, t0, solver->y);
  flag &= CVodeSetUserData(solver->cvode_mem, solver->params);
  flag &= CVDense(solver->cvode_mem, solver->odeModel->N);
  flag &= CVDlsSetDenseJacFn(solver->cvode_mem, solver->odeModel->vf_jac);
  flag &= CVodeSStolerances(solver->cvode_mem, ODE_SOLVER_REL_ERR, ODE_SOLVER_ABS_ERR);
  flag &= CVodeSetMaxNumSteps(solver->cvode_mem, ODE_SOLVER_MX_STEPS);
  if (flag!=CV_SUCCESS) {
    fprintf(stderr,"[CV] ode_solver_init failed, flag=%i\n",flag);
  }
}

void ode_solver_reinit(ode_solver* solver, const double t0,  double* y0, int lenY,  const double* p, int lenP ){
  int i,flag;
  
  /* Get initial conditions */
  if(y0 != 0){
    if( lenY != solver->odeModel->N ){
      fprintf(stderr,"ode_solver_init: lenY must be equal %d, the number of variables in the ode model.\n",solver->odeModel->N);
      return ;
    }
    NV_DATA_S(solver->y) = y0;
  }
  else {
    NV_DATA_S(solver->y) = solver->odeModel->v;
  }
  
  
  /* re-initialise */
  flag = CVodeReInit(solver->cvode_mem, t0, solver->y);
  
  /* Get parameters */
  if (p != 0){
    if (lenP != solver->odeModel->P) {
      fprintf(stderr,"ode_solver_init: lenP must be equal %d, the number of parameters in the ode model.\n",solver->odeModel->P);
      return ;
    }
    int P = solver->odeModel->P;
    for(i = 0; i < P; i++)
      solver->params[i] = p[i];
		
    flag = CVodeSetUserData(solver->cvode_mem, solver->params);
  }
  if (flag!=CV_SUCCESS) {
    fprintf(stderr,"[CV] ode_solver_init failed, flag=%i\n",flag);
  }
}

/* yS0 in row major order */
void ode_solver_init_sens(ode_solver* solver,  double* yS0, int lenP, int lenY){
  /* change: lenP may now be shorter than P if you don't want to
   * calculate all sensitivities.  This is important for the correct
   * treatment of input parameters, i.e. experimental conditions u
   * (since we don't require sensitivities with repsect to them).
   */
  int i,flag;
  
  int N = solver->odeModel->N;
  //int P = solver->odeModel->P;
  
  if (solver->odeModel->vf_sens == 0) {
    fprintf(stderr,"ode_solver_init_sens: no sensitivities defined for this model.\n");
    return;
  }
  
  solver->yS = N_VCloneVectorArrayEmpty_Serial(lenP, solver->y);						/* alloc */
  double tmp[N];
  
  if(yS0 !=0 ){
    if ( (lenY != N) ) {
      fprintf(stderr,"ode_solver_init_sens: lenY must be equal to %d the number of parameters and variables in the ode model.\n",solver->odeModel->N);
      return ;
    }
    
    for(i = 0; i < lenP; i++)
      NV_DATA_S(solver->yS[i]) = &yS0[i*lenY];
  }
  else{
    
    for (i=0; i < N ; i++)
      tmp[i] = 0.0;
    
    for (i = 0; i < lenP; i++)
      NV_DATA_S(solver->yS[i]) = tmp;
  }
  
  flag = CVodeSensInit1(solver->cvode_mem, lenP, CV_STAGGERED1, solver->odeModel->vf_sens, solver->yS);
  flag = CVodeSetSensErrCon(solver->cvode_mem, TRUE);	
  flag = CVodeSensEEtolerances(solver->cvode_mem);
  
  /* set parameters scale for error corection */
  double scale_p[lenP];
  
  for (i=0; i<lenP; i++) {
    /* order of magnitude can be found from the actual parameter value */
    /* if zero then set to the default relative error */
    if (solver->params[i] != 0.0){
      scale_p[i] = solver->params[i];
    }
    else {
      scale_p[i] = ODE_SOLVER_REL_ERR;
    }
    
  }
  flag = CVodeSetSensParams(solver->cvode_mem, solver->params, scale_p, NULL);
  if (flag!=CV_SUCCESS) {
    fprintf(stderr,"[CV] ode_solver_init failed, flag=%i\n",flag);
  }
}

void ode_solver_disable_sens(ode_solver* solver){

  CVodeSensToggleOff(solver->cvode_mem);
  
}

/* yS0 in row major order */
void ode_solver_reinit_sens(ode_solver* solver,  double* yS0, int lenP, int lenY){
  /* change: lenP may now be shorter than P if you dpn't want to
   * calculate all sensitivities.  This is important for the correct
   * treatment of input parameters, i.e. experimental conditions u
   * (since we don't require sensitivities with repsect to them).
   */

  int i,flag;
  int N = solver->odeModel->N;
  //int P = solver->odeModel->P;
  
  if (solver->odeModel->vf_sens == 0) {
    fprintf(stderr,"ode_solver_init_sens: no sensitivities defined for this model.\n");
    return;
  }
  
  double tmp[N];
  if(yS0 !=0 ){
    if ( lenY != N ) {
      fprintf(stderr,"ode_solver_init_sens: lenY must be equal to %d the number of variables in the ode model.\n",solver->odeModel->N);
      return ;
    }
    
    for(i = 0; i < lenP; i++)
      NV_DATA_S(solver->yS[i]) = &yS0[i*lenY];
  }
  else{
    
    for (i=0; i < N ; i++)
      tmp[i] = 0.0;
    
    for (i = 0; i < lenP; i++)
      NV_DATA_S(solver->yS[i]) = tmp;
  }
  
  /* set parameters scale for error corection */
  double scale_p[lenP];
  
  for (i=0; i<lenP; i++) {
    /* order of magnitude can be found from the actual parameter value */
    /* if zero then set to the default relative error */
    if (solver->params[i] != 0.0){
      scale_p[i] = solver->params[i];
    }
    else {
      scale_p[i] = ODE_SOLVER_REL_ERR;
    }
    
  }
  
  flag = CVodeSensReInit(solver->cvode_mem, CV_STAGGERED1, solver->yS);
  flag = CVodeSetSensParams(solver->cvode_mem, solver->params, scale_p, NULL);
  if (flag!=CV_SUCCESS){
    fprintf(stderr,"reinit_sens failed.\n");
    exit(-1);
  }
}

int ode_solver_solve(ode_solver* solver, const double t, double* y, double* tout){
  
  double lTout;
  
  /* Advance the solution */
  NV_DATA_S(solver->y) = y;
  int flag = CVode(solver->cvode_mem,t, solver->y, &lTout, CV_NORMAL);
  
  if(tout)
    *tout = lTout;
  
  return flag;
}

/* yS in row major order */
void ode_solver_get_sens(ode_solver* solver, double t, double* yS){
  int i, N_s, P_s;
  
  N_s = solver->odeModel->N;
  P_s = solver->odeModel->P;
  
  /* return the sensitivities to the return parameter yS */
  if ( solver->odeModel->vf_sens != 0 ){
    for(i = 0; i < P_s; i++)
      NV_DATA_S(solver->yS[i]) = &yS[i*N_s];
    
    CVodeGetSens(solver->cvode_mem, &t, solver->yS);
  }
}

void ode_solver_get_func(ode_solver* solver, const double t, double* y, double* fy){
  
  NV_DATA_S(solver->y) = y;
  
  /* evaluate and return the functions to the return parameter fy */
  solver->odeModel->vf_func(t, solver->y, fy ,(void *) solver->params);
  
}

/* yS and fyS in row major order */
void ode_solver_get_func_sens(ode_solver* solver, const double t,  double* y,  double* yS, double* fyS){
  int i;
  int N = solver->odeModel->N;
  int P = solver->odeModel->P;
  
  NV_DATA_S(solver->y) = y;
  
  if ( solver->odeModel->vf_sens != 0 ){
    for(i = 0; i < P; i++)
      NV_DATA_S(solver->yS[i]) = &yS[N*i];
    
    /* evaluate and return the functions sensitivities to the return parameter fyS */
    if ( solver->odeModel->vf_func_sens != 0  )
      solver->odeModel->vf_func_sens(t, solver->y, solver->yS , fyS, (void *) solver->params);
  }
}

void ode_solver_print_stats(const ode_solver* solver, FILE* outF){
  long int nst;
  long int nfe, nsetups, nni, ncfn, netf;
  long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
  long int nje, nfeLS;
  int flag;
  
  void* cvode_mem = solver->cvode_mem;
  
  flag = CVodeGetNumSteps(cvode_mem, &nst);
  
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  
  
  if (solver->yS != 0) {
    flag = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
	
    flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
    
    flag = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
    
    flag = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
    
    flag = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
    
    flag = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
    
  }
  
  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  
  
  fprintf(outF,"\n# Solver Statistics\n\n");
  fprintf(outF,"# Steps            = %5ld\n\n", nst);
  fprintf(outF,"# RhsEvals         = %5ld\n",   nfe);
  fprintf(outF,"# ErrTestFails     = %5ld   LinSolvSetups        = %5ld\n", netf, nsetups);
  fprintf(outF,"# NonlinSolvIters  = %5ld   NonlinSolvConvFails  = %5ld\n", nni, ncfn);
  
  if(solver->yS != 0) {
    fprintf(outF,"\n# Sensitivities Statistics\n");
    fprintf(outF,"# SensRhsEvals     = %5ld   RhsEvals             = %5ld\n", nfSe, nfeS);
    fprintf(outF,"# ErrTestFails     = %5ld   LinSolvSetups        = %5ld\n", netfS, nsetupsS);
    fprintf(outF,"# NonlinSolvIters  = %5ld   NonlinSolvConvFails  = %5ld\n", nniS, ncfnS);
  }
  
  fprintf(outF,"\n# Jacobian Statistics\n");
  fprintf(outF,"# JacEvals  = %5ld    RhsEvals  = %5ld\n", nje, nfeLS);
  if (flag!=CV_SUCCESS) {
    fprintf(stderr,"[CV] ode_solver_init failed, flag=%i\n",flag);
  }
}


void ode_solver_free(ode_solver* solver){
  free(solver->params);
  CVodeFree(&(solver->cvode_mem));
  if (solver->yS != 0) {
    N_VDestroyVectorArray_Serial(solver->yS, solver->odeModel->P);
  }
  N_VDestroy_Serial(solver->y);
  free(solver);
}
