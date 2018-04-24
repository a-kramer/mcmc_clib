#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "../ode/ode_model.h"
#include "../mcmc/model_parameters_smmala.h"
#include "read_data_hdf5.h"

herr_t load_stdv_block(hid_t g_id, const char *name, const H5L_info_t *info, void *op_data){
  ode_model_parameters *mp=op_data;
  int rank;
  herr_t status=H5LTget_dataset_ndims(g_id, name, &rank);
  hsize_t *size;
  static int index=0;
  size=malloc(sizeof(size_t)*rank);
  status&=H5LTget_dataset_info(g_id,name,size,NULL,NULL);
  assert(status>=0); // I think that hdf5 functions return negative error codes
  gsl_matrix *stdv_block;
  stdv_block=gsl_matrix_alloc(size[0],size[1]);
  status&=H5LTread_dataset_double(g_id, name, stdv_block->data);
  if (index < mp->size->C){
    mp->E[index]->sd_data_block=stdv_block;
    index++;
  }else{
    fprintf(stderr,"[load_stdv_block] error: experiment index is larger than number ofexperiments (simulation units).\n");
    exit(-1);
  }
  return status>=0?0:-1;  
}

herr_t load_data_block(hid_t g_id, const char *name, const H5L_info_t *info, void *op_data){
  ode_model_parameters *mp=op_data;
  int rank;
  herr_t status=H5LTget_dataset_ndims(g_id, name, &rank);
  hsize_t *size;
  size=malloc(sizeof(size_t)*rank);
  status&=H5LTget_dataset_info(g_id,name,size,NULL,NULL);
  assert(status>=0); // I think that hdf5 functions return negative error codes
  gsl_matrix *data_block;
  data_block=gsl_matrix_alloc(size[0],size[1]);
  status&=H5LTread_dataset_double(g_id, name, data_block->data);
  int index,major,minor;
  status&=H5LTget_attribute_int(g_id,name,"index",&index);
  status&=H5LTget_attribute_int(g_id,name,"major",&major);
  status&=H5LTget_attribute_int(g_id,name,"minor",&minor);
  hsize_t T,U,D; // measurement time, lenght of input vector, number of unknown parameters
  status&=H5LTget_attribute_info(g_id,name,"time",&T,NULL,NULL);
  printf("[load_experiment] loading simulation unit %i; Experiment %i.%i with %lli measurements\n",index,major,minor,T);
  gsl_vector *time, *input;
  time=gsl_vector_alloc(T);
  H5LTget_attribute_double(g_id,name,"time",time->data);
  H5LTget_attribute_info(g_id,name,"input",&U,NULL,NULL);
  input=gsl_vector_alloc(U);
  H5LTget_attribute_double(g_id,name,"input",input->data);
  // normalisation properties
  int NormaliseByExperiment;
  int NormaliseByTimePoint;
  gsl_vector_int *NormaliseByOutput;
  herr_t attr_err;
  attr_err=H5LTget_attribute_int(g_id, name,"NormaliseByExperiment",&NormaliseByExperiment,1);
  if (attr_err<0) NormaliseByExperiment=-1;
  attr_err=H5LTget_attribute_int(g_id, name,"NormaliseByTimePoint",&NormaliseByTimePoint,1);
  if (attr_err<0) NormaliseByTimePoint=-1;
  hsize_t nO;  
  attr_err=H5LTget_attribute_info(g_id,name,"NormaliseByOutput",&nO,NULL,NULL);
  if (attr_err<0) {
    NormaliseByOutput=NULL;
  } else {
    NormaliseByOutput=gsl_vector_int_alloc(nO);
    attr_err=H5LTget_attribute_int(g_id,name,"NormaliseByOutput",NormaliseByOutput->data,nO);
    assert(attr_err>=0);
  }
  
  if (index < mp->size->C){
    mp->E[index]->data_block=data_block;
    mp->E[index]->input_u=input;
    mp->E[index]->t=time;
    mp->E[index]->input_u=input;
    mp->E[index]->NormaliseByTimePoint=NormaliseByTimePoint;
    mp->E[index]->NormaliseByTimePoint=NormaliseByTimePoint;
    mp->E[index]->NormaliseByTimePoint=NormaliseByTimePoint;
    mp->size->T=GSL_MAX(mp->size->T,mp->E[index]->t->size);
  }else{
    fprintf(stderr,"[load_experiment_block] error: experiment index is larger than number ofexperiments (simulation units).\n");
    exit(-1);
  }

  return status>=0?0:-1;
}

/* this function loads the experimental data from an hdf5 file the
 * second argument is called model_parameters because it stores
 * everything that must be part of the "ODE model" for a specific mcmc
 * algorithm to run. It contains the data but also allocated space for
 * simulation results and much more
 */
int read_data(const char *file, void *model_parameters){
  ode_model_parameters *mp=model_parameters;
  ode_solver *solver=mp->solver;
  ode_model *M=solver->odeModel;
  // these values should alrfeady have been set
  //mp->size->N=ode_model_getN(M);
  //mp->size->P=ode_model_getP(M);
  //mp->size->F=ode_model_getF(M);
  
  herr_t status;
  hid_t file_id = H5Fopen(file,H5F_ACC_RDONLY,H5P_DEFAULT);
  hid_t data_group_id = H5Gopen2(file_id,"/data",H5P_DEFAULT);
  hid_t stdv_group_id = H5Gopen2(file_id,"/sd_data",H5P_DEFAULT);
  assert(data_group_id>0 && stdv_group_id>0);
  mp->size->C;
  hsize_t idx,nE;
  status&=H5Gget_num_objs(data_group_id, &nE); //number of experiments
  mp->size->C=(int) nE;
  mp->size->T=0;
  init_E(mp);
  idx=0;
  status=H5Literate(data_group_id, H5_INDEX_NAME, H5_ITER_INC, &idx, load_data_block, mp);
  if (status<0) fprintf(stderr,"[read_data] iteration over data group was not successful.\n");
  idx=0;
  status=H5Literate(stdv_group_id, H5_INDEX_NAME, H5_ITER_INC, &idx, load_data_block, mp);
  if (status<0) fprintf(stderr,"[read_data] iteration over stdv group was not successful.\n");  
  H5Gclose(data_group_id);
  H5Gclose(stdv_group_id);
  hid_t prior_group_id=H5Gopen2(file_id,"/prior",H5P_DEFAULT);
  int D;
  status&=H5LTget_dataset_info(prior_group_id,"mu",&D,NULL,NULL);

  mp->prior_mu=gsl_vector_alloc(D);
  status&=H5LTread_dataset_double(prior_group_id, "mu", mp->prior_mu->data);
  
  mp->Sigma=gsl_matrix_alloc(D,D);
  status&=H5LTread_dataset_double(prior_group_id, "Sigma", mp->prior_mu->data);

  H5Gclose(prior_group_id);
  H5Fclose(file_id);
  // determine sizes:
  mp->size->U=mp->E[0]->input_u->size;
  mp->size->D=D;
  ode_model_parameters_link(mp);
  return (int) status;
}
