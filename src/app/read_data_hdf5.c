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
#include <gsl/gsl_permutation.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "../ode/ode_model.h"
#include "../mcmc/model_parameters_smmala.h"
#include "read_data_hdf5.h"
#include "normalisation_sd.h"


void* h5_to_gsl_int(hid_t g_id, const char *obj_name, const char* attr_name){
  int rank;
  hsize_t size[2];
  herr_t status=0;
  gsl_vector_int *v;
  gsl_matrix_int *m;
  /*attributes always have rank==1, at least when the H5LT API is used
    to make them*/
  if (attr_name){
    status|=H5LTget_attribute_info(g_id,obj_name,attr_name,size,NULL,NULL);
    v=gsl_vector_int_alloc(size[0]);
    status|=H5LTget_attribute_int(g_id,obj_name,attr_name,v->data);
  } else {
    H5LTget_dataset_ndims(g_id, obj_name, &rank);
    H5LTget_dataset_info(g_id, obj_name, size, NULL, NULL);
    if (rank==1){
      v=gsl_vector_int_alloc(size[0]);
      H5LTread_dataset_int(g_id, obj_name, v->data);
      return v;
    } else if (rank==2){
      m=gsl_matrix_int_alloc(size[0],size[1]);
      H5LTread_dataset_int(g_id, obj_name, m->data);
      return m;
    } else {
      //error
      fprintf(stderr,"[%s] hdf5 data object has an invalid number of dimensions: %i.\n",__func__,rank);
      return NULL;
    }
  }
  return v;
}


void* h5_to_gsl(hid_t g_id, const char *obj_name, const char* attr_name){
  int rank;
  hsize_t size[2];
  herr_t status=0;
  gsl_vector *v;
  gsl_matrix *m;
  /*attributes always have rank==1, at least when the H5LT API is used
    to make them*/
  if (attr_name){
    status|=H5LTget_attribute_info(g_id,obj_name,attr_name,size,NULL,NULL);
    v=gsl_vector_alloc(size[0]);
    status|=H5LTget_attribute_double(g_id,obj_name,attr_name,v->data);
  } else {
    H5LTget_dataset_ndims(g_id, obj_name, &rank);
    H5LTget_dataset_info(g_id, obj_name, size, NULL, NULL);
    if (rank==1){
      v=gsl_vector_alloc(size[0]);
      H5LTread_dataset_double(g_id, obj_name, v->data);
      return v;
    } else if (rank==2){
      m=gsl_matrix_alloc(size[0],size[1]);
      H5LTread_dataset_double(g_id, obj_name, m->data);
      return m;
    } else {
      //error
      fprintf(stderr,"[%s] hdf5 data object has an invalid number of dimensions: %i.\n",__func__,rank);
      return NULL;
    }
  }
  return v;
}

/* loads each standard deviation dataset from the hdf5 data file
 * and saves the information in the experiment structure indexed as
 * annotated in the hdf5 dataset. Both the data and stdv block carry
 * an index attribute and are then saved in @code E[index]@. 
 */
herr_t /*hdf5 error type, undocumented*/
load_stdv_block(hid_t g_id, /*group id (hdf5)*/
		const char *name, /*dataset name*/
		const H5L_info_t *info, /*info struct, as defined in hdf5 API*/
		void *op_data)/*model parameter struct, contains experiment array*/{
  ode_model_parameters *mp=op_data;
  int rank;
  herr_t status=H5LTget_dataset_ndims(g_id, name, &rank);
  hsize_t *size;
  int index=-1;
  size=malloc(sizeof(size_t)*rank);
  status|=H5LTget_dataset_info(g_id,name,size,NULL,NULL);
  assert(status>=0); // I think that hdf5 functions return negative error codes
  gsl_matrix *stdv_block;
  stdv_block=gsl_matrix_alloc(size[0],size[1]);
  status|=H5LTread_dataset_double(g_id, name, stdv_block->data);

  int major, minor;
  //printf("[load_data_block] checking experiment's index: ");
  status|=H5LTget_attribute_int(g_id,name,"index",&index); //printf(" %i ",index); assert(status>=0);
  status|=H5LTget_attribute_int(g_id,name,"major",&major); //printf("%i.",major); assert(status>=0);
  status|=H5LTget_attribute_int(g_id,name,"minor",&minor); //printf("%i",minor);
  assert(status>=0);
  assert(index>=0);
  
  /* This is a debugging comparison by name using scanf. Not necessary.
   */
  //int i,nout=0;
  //printf("[load_stdv_block] the current hdf5 «name» is: %s\n",name);
  //nout=sscanf(name,"sd_data_block_%d",&i);
  //printf("[load_stdv_block] nout=%d; the retrieved stdv index is: %i\n",nout,index);
  //assert(nout==1);
  //assert(i==index);
  
  if (index < mp->size->C && index >= 0){
    assert(mp->E[index]->sd_data_block==NULL);
    mp->E[index]->sd_data_block=stdv_block;
  }else{
    fprintf(stderr,"[load_stdv_block] error: experiment index is larger than number ofexperiments (simulation units).\n");
    exit(-1);
  }
  free(size);
  return status>=0?0:-1;  
}

/* loads each experimental data block from the hdf5 file and saves the
 * information in the experiment structure indexed as annotated. Both
 * the data and stdv block carry an index attribute and are then saved
 * in @code E[index]@.  
 */
herr_t /*hdf5 error type, undocumented*/
load_data_block(hid_t g_id,
		const char *name,
		const H5L_info_t *info,
		void *op_data)/*ode model parameters, contains experiment array E[]*/{
  ode_model_parameters *mp=op_data;
  int rank;
  herr_t status=H5LTget_dataset_ndims(g_id, name, &rank);
  hsize_t *size;
  size=malloc(sizeof(hsize_t)*rank);
  status|=H5LTget_dataset_info(g_id,name,size,NULL,NULL);
  assert(status>=0); // I think that hdf5 functions return negative error codes, but it's undocumented
  gsl_matrix *data_block;
  //printf("[load_data_block] loading data of size %lli×%lli.\n",size[0],size[1]);
  fflush(stdout);
  data_block=gsl_matrix_alloc((size_t) size[0],(size_t) size[1]);
  assert(data_block!=NULL && data_block->data !=NULL);  
  status|=H5LTread_dataset_double(g_id, name, data_block->data);
  assert(status>=0);
  //printf("[load_data_block] gsl_matrix(%li,%li) loaded.\n",data_block->size1,data_block->size2);
  fflush(stdout);
  int index, major, minor, lflag;
  //printf("[load_data_block] checking experiment's index: ");
  status|=H5LTget_attribute_int(g_id,name,"index",&index); //printf(" %i ",index); assert(status>=0);
  status|=H5LTget_attribute_int(g_id,name,"major",&major); //printf("%i.",major); assert(status>=0);
  status|=H5LTget_attribute_int(g_id,name,"minor",&minor); //printf("%i",minor); assert(status>=0);
  status|=H5LTget_attribute_int(g_id,name,"LikelihoodFlag",&lflag); //printf(" %i ",lflag);
  //printf("\n");
  assert(status>=0);
  hsize_t T,U; // measurement time, lenght of input vector, number of unknown parameters
  status|=H5LTget_attribute_ndims(g_id,name,"time",&rank);
  //printf("[load_data_block] checking experiment's measurement times(%i).\n",rank);
  fflush(stdout);
  assert(status>=0);
  H5T_class_t type_class;
  size_t      type_size;
  gsl_vector *time, *input;
  
  //  if (rank==1){    
  status|=H5LTget_attribute_info(g_id,name,"time",&T,&type_class,&type_size);
  assert(status>=0);    
  //printf("[load_data_block] loading simulation unit %i; Experiment %i.%i with %lli measurements\n",index,major,minor,T);    
  fflush(stdout);
  time=gsl_vector_alloc((size_t) T);
  H5LTget_attribute_double(g_id,name,"time",time->data);
  //}
  H5LTget_attribute_info(g_id,name,"input",&U,&type_class,&type_size);
  //printf("[load_data_block] input size: %i.\n",(int) U);
  input=gsl_vector_alloc((size_t) U);
  H5LTget_attribute_double(g_id,name,"input",input->data);
  
  // normalisation properties
  int NormaliseByExperiment;
  int NormaliseByTimePoint;
  gsl_vector_int *NormaliseByOutput=NULL;
  herr_t attr_err;
  hid_t d_id=H5Dopen2(g_id, name, H5P_DEFAULT);
  if (H5LTfind_attribute(d_id,"NormaliseByExperiment")){
    attr_err=H5LTget_attribute_int(g_id, name,"NormaliseByExperiment",&NormaliseByExperiment);
    assert(attr_err>=0);
  } else {
    NormaliseByExperiment=-1;
  }  
  if (H5LTfind_attribute(d_id,"NormaliseByTimePoint")){
    attr_err=H5LTget_attribute_int(g_id, name,"NormaliseByTimePoint",&NormaliseByTimePoint);
    assert(attr_err>=0);
  } else{
    NormaliseByTimePoint=-1;
  }
  hsize_t nO;
  if (H5LTfind_attribute(d_id,"NormaliseByOutput")){
    attr_err=H5LTget_attribute_info(g_id,name,"NormaliseByOutput",&nO,&type_class,&type_size);    
    if (attr_err>=0) {
      NormaliseByOutput=gsl_vector_int_alloc(nO);
      attr_err=H5LTget_attribute_int(g_id,name,"NormaliseByOutput",NormaliseByOutput->data);
      assert(attr_err>=0);
    } else {
      fprintf(stderr,"[load_data_block] read error for NormaliseByOutput vector.\n");
    }
  } else {
    NormaliseByOutput=NULL;
  }
  
  if (index < mp->size->C && index >= 0){
    mp->E[index]->index[0]=index;
    mp->E[index]->index[1]=major;
    mp->E[index]->index[2]=minor;    
    mp->E[index]->lflag=lflag;
    assert(mp->E[index]->data_block==NULL);
    mp->E[index]->data_block=data_block;
    mp->E[index]->input_u=input;
    mp->E[index]->t=time;
    mp->E[index]->input_u=input;
    mp->E[index]->NormaliseByExperiment=NormaliseByExperiment;
    mp->E[index]->NormaliseByTimePoint=NormaliseByTimePoint;
    mp->E[index]->NormaliseByOutput=NormaliseByOutput;
    mp->size->T=GSL_MAX(mp->size->T,mp->E[index]->t->size);
  }else{
    fprintf(stderr,"[load_experiment_block] error: experiment index is larger than number ofexperiments (simulation units).\n");
    exit(-1);
  }
  free(size);
  return status>=0?0:-1;
}


/* H5Literate function to load prior parameters @code op_data@ is the model
 * parameters struct. The function copies mu and one of: sigma, Sigma, Precision into the struct
 * Checks what type of prior has been set.
 */
int /*returns GSL error codes*/
load_prior(hid_t g_id, void *op_data) /*ode model parameters, contains prior parameter struct*/{
  ode_model_parameters *mp=op_data;
  hsize_t D;
  herr_t status=0;
  int gsl_err=GSL_SUCCESS;
  status=H5LTget_dataset_info(g_id,"mu",&D,NULL,NULL);
  //printf("[read_data] MCMC sampling will be performed on the first %lli parameters.\n",D);
  mp->size->D=(int) D;
  mp->prior->mu=gsl_vector_alloc((size_t) D);
  status|=H5LTread_dataset_double(g_id, "mu", mp->prior->mu->data);
  
  if (H5LTfind_dataset(g_id,"Sigma")){
    mp->prior->Sigma_LU=gsl_matrix_alloc(D,D);
    status|=H5LTread_dataset_double(g_id, "Sigma", mp->prior->Sigma_LU->data);
    mp->prior->type=(PRIOR_IS_GAUSSIAN | PRIOR_IS_MULTIVARIATE | PRIOR_SIGMA_GIVEN);
    mp->prior->p=gsl_permutation_alloc((size_t) D);
    gsl_err|=gsl_linalg_LU_decomp(mp->prior->Sigma_LU, mp->prior->p, &(mp->prior->signum));
  } else if (H5LTfind_dataset(g_id,"sigma")){
    mp->prior->sigma=gsl_vector_alloc(D);
    status|=H5LTread_dataset_double(g_id, "sigma", mp->prior->sigma->data);
    mp->prior->type=(PRIOR_IS_GAUSSIAN);
  } else if (H5LTfind_dataset(g_id,"Precision")){
    mp->prior->inv_cov=gsl_matrix_alloc(D,D);
    status|=H5LTread_dataset_double(g_id, "Precision", mp->prior->inv_cov->data);
    mp->prior->type=(PRIOR_IS_GAUSSIAN | PRIOR_IS_MULTIVARIATE | PRIOR_PRECISION_GIVEN);
  } else {
    fprintf(stderr,"[read_data] not enough prior information given. (Specify either Sigma, sigma, or inverse_covariance)\n");
    exit(-1);
  }
  //  assert(gsl_err==GSL_SUCCESS);
  return gsl_err;
}


/* loads events if present
*/
herr_t /*hdf5 error type, undocumented*/
load_event_block(hid_t g_id, /*group id (hdf5)*/
		const char *name, /*dataset name*/
		const H5L_info_t *info, /*info struct, as defined in hdf5 API*/
		void *op_data)/*model parameter struct, contains experiment array*/{
  ode_model_parameters *mp=op_data;
  int rank;
  herr_t status=H5LTget_dataset_ndims(g_id, name, &rank);
  hsize_t size[2]; // rank is always 1 or 2 in all of our cases;
  int index=-1;
  //size=malloc(sizeof(size_t)*rank);
  status|=H5LTget_dataset_info(g_id,name,size,NULL,NULL);
  assert(status>=0); // I think that hdf5 functions return negative error codes
  gsl_matrix *event_value=h5_to_gsl(g_id,name,NULL);
  gsl_vector_int *event_target=h5_to_gsl_int(g_id,name,"Target");
  gsl_vector *event_time=h5_to_gsl(g_id,name,"Time");
  gsl_vector_int *event_type=h5_to_gsl(g_id,name,"Type");
  
  int major, minor;
  status|=H5LTget_attribute_int(g_id,name,"index",&index);
  status|=H5LTget_attribute_int(g_id,name,"major",&major);
  status|=H5LTget_attribute_int(g_id,name,"minor",&minor);
  assert(status>=0);
  assert(index>=0);
  experiment *E;
  event_t *event;
  if (index < mp->size->C && index >= 0){
    E=mp->E[index];
    event=add_event(E->event,E->t,event_time,event_type,event_target,event_value);
    assert(E->single!=NULL);
    /* insert events into global time line */
    insert_events(E->t, E->single, event); 
  }else{
    fprintf(stderr,"[%s] error: experiment index is larger than number ofexperiments (simulation units).\n",__func__);
    exit(-1);
  }
  //free(size);
  return status>=0?0:-1;  
}




/* loads the experimental data from an hdf5 file the
 * second argument is called model_parameters because it stores
 * everything that must be part of the "ODE model" for a specific mcmc
 * algorithm to run. It contains the data but also allocated space for
 * simulation results and much more.
 */
int /*casts hdf5 error codes to int*/
read_data(const char *file,
	  void *model_parameters)/*a fairly big struct that contains
	 data, prior parameters, and pre allocated space for
	 simulation results and normalisation calculations*/{
  ode_model_parameters *mp=model_parameters;
  herr_t status;
  hid_t file_id = H5Fopen(file,H5F_ACC_RDONLY,H5P_DEFAULT);
  //printf("[read_data] HDF5 file id=%li.\n",file_id); fflush(stdout);
  hid_t data_group_id = H5Gopen2(file_id,"/data",H5P_DEFAULT);
  //printf("[read_data] HDF5 data group id=%li.\n",data_group_id); fflush(stdout);
  hid_t stdv_group_id = H5Gopen2(file_id,"/stdv",H5P_DEFAULT);
  //printf("[read_data] HDF5 standard deviation group id=%li.\n",stdv_group_id); fflush(stdout);
  assert(file_id>0 && data_group_id>0 && stdv_group_id>0);
  hsize_t idx,nE;
  int gsl_err;
  status|=H5Gget_num_objs(data_group_id, &nE); //number of experiments
  mp->size->C=(int) nE;
  mp->size->T=0;
  init_E(mp);
  printf("[read_data] experiments initialised.\n"); fflush(stdout);
  idx=0;
  status=H5Literate(data_group_id, H5_INDEX_NAME, H5_ITER_INC, &idx, load_data_block, mp);
  if (status<0) {
    fprintf(stderr,"[%s] iteration over data group was not successful.\n",__func__);
    fflush(stderr);
  } else {
    printf("[%s] data group read.\n",__func__);
  }
  fflush(stdout);
  idx=0;
  status=H5Literate(stdv_group_id, H5_INDEX_NAME, H5_ITER_INC, &idx, load_stdv_block, mp);
  if (status<0){
    fprintf(stderr,"[%s] iteration over stdv group was not successful.\n",__func__);
  }else{
    printf("[%s] standard deviations read: «stdv» group.\n",__func__);
  }
  H5Gclose(data_group_id);
  H5Gclose(stdv_group_id);
  size_t i;
  /* read event tables */
  if (H5Lexists(file_id,"/event",H5P_DEFAULT)>0){
    hid_t event_group_id = H5Gopen2(file_id,"/event",H5P_DEFAULT);
    idx=0;
    status = H5Literate(event_group_id, H5_INDEX_NAME, H5_ITER_INC, &idx, load_event_block, mp);
    assert(status>=0);
    for (i=0;i<mp->size->C;i++){
      if (mp->E[i]->event)
	mp->E[i]->before_t=convert_to_array(mp->E[i]->t->size,mp->E[i]->single);
    }
  }
  hid_t prior_group_id=H5Gopen2(file_id,"/prior",H5P_DEFAULT);
  gsl_err = load_prior(prior_group_id,model_parameters);
  assert(gsl_err==GSL_SUCCESS);
  H5Gclose(prior_group_id);
  H5Fclose(file_id);
  // determine sizes:
  mp->size->U=(int) (mp->E[0]->input_u->size);
  //printf("[read_data] Simulations require %i known input parameters.\n",mp->size->U);
  return (int) status;
}
