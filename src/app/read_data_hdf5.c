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
#include "../sb/h5block.h"
#include "../ode/ode_model.h"
#include "../mcmc/model_parameters_smmala.h"
#include "read_data_hdf5.h"
#include "normalisation_sd.h"

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
  herr_t status=0;
  int index=-1;
  gsl_matrix *stdv_block=h5_to_gsl(g_id,name,NULL);
  int major, minor;
  status|=H5LTget_attribute_int(g_id,name,"index",&index); //printf(" %i ",index); assert(status>=0);
  status|=H5LTget_attribute_int(g_id,name,"major",&major); //printf("%i.",major); assert(status>=0);
  status|=H5LTget_attribute_int(g_id,name,"minor",&minor); //printf("%i",minor);
  assert(status>=0);
  assert(index>=0);
  if (index < mp->size->C && index >= 0){
    assert(mp->E[index]->sd_data_block==NULL);
    mp->E[index]->sd_data_block=stdv_block;
  }else{
    fprintf(stderr,"[load_stdv_block] error: experiment index is larger than number ofexperiments (simulation units).\n");
    exit(-1);
  }
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
  //int rank;
  herr_t status=0;
  gsl_matrix *data_block=h5_to_gsl(g_id,name,NULL);
  assert(data_block!=NULL && data_block->data !=NULL);  
  //printf("[load_data_block] gsl_matrix(%li,%li) loaded.\n",data_block->size1,data_block->size2);
  //  fflush(stdout);
  int index, major, minor, lflag;
  status|=H5LTget_attribute_int(g_id,name,"index",&index);
  status|=H5LTget_attribute_int(g_id,name,"major",&major);
  status|=H5LTget_attribute_int(g_id,name,"minor",&minor);
  status|=H5LTget_attribute_int(g_id,name,"LikelihoodFlag",&lflag);
  assert(status>=0);
  gsl_vector *time=h5_to_gsl(g_id,name,"time");
  gsl_vector *input=h5_to_gsl(g_id,name,"input");
  gsl_vector *y0=h5_to_gsl(g_id,name,"InitilValue");
  // normalisation properties
  int NormaliseByExperiment;
  int NormaliseByTimePoint;
  gsl_vector_int *NormaliseByOutput=NULL;
  herr_t attr_err;
  H5T_class_t type_class;
  size_t type_size;

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
    mp->E[index]->init_y=y0;
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


/* loads events 
*/
herr_t /*hdf5 error type, undocumented*/
load_event_block
(hid_t g_id, /*group id (hdf5)*/
 const char *name, /*dataset name*/
 const H5L_info_t *info, /*info struct, as defined in hdf5 API*/
 void *op_data)/*model parameter struct, contains experiment array*/
{
  ode_model_parameters *mp=op_data;
  assert(mp->model_p.name);
  assert(mp->model_x.name);

  int rank;
  // printf("[%s] getting dataset size («%s»).\n",__func__,name); fflush(stdout);
  herr_t status=H5LTget_dataset_ndims(g_id, name, &rank);
  assert(rank<=2 && rank>0);
  hsize_t size[2]; // rank is always 1 or 2 in all of our cases, so 2 always works;
  status|=H5LTget_dataset_info(g_id,name,size,NULL,NULL);
  // printf("[%s] size of «%s» is %lli × %lli.\n",__func__,name,size[0],size[1]); fflush(stdout);
  assert(status>=0); // I think that hdf5 functions return negative error codes
  
  //char *ExperimentName=h5_to_char(g_id, name, "ExperimentName");
  // printf("[%s] affected experiment «%s».\n",__func__,ExperimentName);
  //fflush(stdout);
  char *TargetName=h5_to_char(g_id, name, "TargetName");
  printf("[%s] affected variable «%s».\n",__func__,TargetName);
  fflush(stdout);
  gsl_matrix *value=h5_to_gsl(g_id,name,NULL);
  assert(value);
  gsl_vector *time=h5_to_gsl(g_id,name,"Time");
  assert(time);
  gsl_vector_fprintf(stdout,time,"%g"); fflush(stdout);
  gsl_vector_int *type=h5_to_gsl_int(g_id,name,"op");
  assert(type);
  gsl_vector_int *effect=h5_to_gsl_int(g_id,name,"Effect");
  assert(effect);  
  gsl_vector_int *target=event_find_targets
    ((effect_t*) effect->data,
     TargetName, effect->size,
     mp->model_p.name,mp->model_p.size,
     mp->model_x.name,mp->model_x.size
     );
  assert(target);
  //printf("[%s] event targets:\n",__func__);
  gsl_vector_int_fprintf(stdout,target,"%i");
  int index=-1;
  int major, minor;
  status|=H5LTget_attribute_int(g_id,name,"index",&index);
  status|=H5LTget_attribute_int(g_id,name,"AffectsMajorIndex",&major);
  status|=H5LTget_attribute_int(g_id,name,"AffectsMinorIndex",&minor);
  //printf("[%s] Experiment: %i (%i.%i).\n",__func__,index,major,minor);
  fflush(stdout);
  assert(status>=0);
  assert(index>=0);
  experiment *E;
  event_t *event;
  if (index < mp->size->C && index >= 0){
    assert(mp->E[index]);
    E=mp->E[index];
    mp->E[index]->event_list=event_list_alloc(3);
    event=event_append(E->event_list,E->t,time,type,effect,target,value);
    /* insert events into global time line E->t */
    event_push(E->single,E->t,event);
  }else{
    fprintf(stderr,"[%s] error: experiment index out of bounds 0<=%i<%i.\n",__func__,index,mp->size->C);
    abort();
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
read_data
(const char *file,
 void *model_parameters)/*a fairly big struct that contains
			  data, prior parameters, and pre allocated space for
			  simulation results and normalisation calculations*/
{
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
  status=H5Gget_num_objs(data_group_id, &nE); //number of experiments
  mp->size->C=(int) nE;
  mp->size->T=0;
  init_E(mp);
  printf("[%s] experiments initialised.\n",__func__); fflush(stdout);
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
    printf("[%s] loading events from file.\n",__func__); fflush(stdout);
    hid_t event_group_id = H5Gopen2(file_id,"/event",H5P_DEFAULT);
    idx=0;
    for (i=0;i<mp->size->C;i++) mp->E[i]->single=calloc(mp->E[i]->t->size,sizeof(event_row_t*));
    status = H5Literate(event_group_id, H5_INDEX_NAME, H5_ITER_INC, &idx, load_event_block, mp);
    printf("[%s] events read.\n",__func__);
    assert(status>=0);
    for (i=0;i<mp->size->C;i++){
      if (mp->E[i]->event_list)
	mp->E[i]->before_t=event_convert_to_arrays(mp->E[i]->t->size,mp->E[i]->single);
      events_are_time_ordered(mp->E[i]->t->size,mp->E[i]->before_t);
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
