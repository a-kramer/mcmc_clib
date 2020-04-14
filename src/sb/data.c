#include "data.h"
#include "re.h"
#include <math.h>

data_t* data_block_alloc(char *BaseName, int major, int minor, gsl_matrix *data, gsl_matrix *stdv, gsl_vector *input, gsl_vector *time, int likelihood_flag){
  assert(time->size == data->size1);
  assert(data->size1 == stdv->size1 && data->size2 == stdv->size2);
  data_t *D=malloc(sizeof(data_t));
  GString *Name=g_string_new(NULL);
  g_string_printf(Name,"E-%i.%i-%s",major,minor,BaseName);
  D->MajorIndex=major;
  D->MinorIndex=minor;
  D->measurement=data;
  D->noise=stdv;
  D->input=input;
  D->time=time;
  D->lflag=likelihood_flag;
  D->Name=Name;
  return D;
}

GPtrArray* /* Data blocks. Each can be simulated with exactly one input vector. */
unwrap_data
(const GPtrArray *DataTable, /* Data tables exactly as they appear in the spreadsheet*/
 const int *lflag, /* the measurement contributes to Likelihood (or not)*/
 const experiment_type *ExperimentType, /* Type: _Dose Response_ or _Time Series_ */
 const GPtrArray *ExperimentName,
 const GPtrArray *input_override, /* default input for each experiment */
 const gsl_vector *default_time, /* default measurement time */
 const sbtab_t *Input,
 const sbtab_t *Output,
 map_t *IdxMap)
{
  int i,j;
  int L;
  gsl_matrix *Y_DATA, *Y_STDV;
  gsl_vector *u_j, *time;
  int nE=DataTable->len;
  sbtab_t *DataTable_j;
  GPtrArray *Data=g_ptr_array_sized_new(nE);
  data_t *D;

  GPtrArray *uid=sbtab_find_column(Input,"!ID",NULL);
  GPtrArray *oid=sbtab_find_column(Output,"!ID",NULL);
  GPtrArray *ErrorName=sbtab_find_column(Output,"!ErrorName",NULL);
  gsl_matrix *input_block;
  gsl_vector_view input_row, time_view;
  gsl_matrix_view data_sub, stdv_sub;
  int nU=uid->len;
  char *E_j;
  assert(lflag);
  for (j=0;j<nE;j++) { // j is the major experiment index
    DataTable_j=g_ptr_array_index(DataTable,j);
    L=table_length(DataTable_j);
    u_j=g_ptr_array_index(input_override,j);
    E_j=g_ptr_array_index(ExperimentName,j);
    if (DataTable_j){
      Y_DATA=sbtab_columns_to_gsl_matrix(DataTable_j,oid,">",1.0);
      Y_STDV=sbtab_columns_to_gsl_matrix(DataTable_j,ErrorName,NULL,INFINITY);      
      switch(ExperimentType[j]){
      case time_series:
	time=sbtab_column_to_gsl_vector(DataTable_j,"!Time");
	map_index(IdxMap,j,0);
	D=data_block_alloc(E_j,j,0,Y_DATA,Y_STDV,u_j,time,lflag[j]);
	g_ptr_array_add(Data,D);
	break;
      case dose_response:
	input_block=gsl_matrix_alloc(L,nU);
	for (i=0;i<L;i++) {
	  input_row=gsl_matrix_row(input_block,i);
	  gsl_vector_memcpy(&(input_row.vector),u_j);
	}
	sbtab_update_gsl_matrix(input_block,DataTable_j,uid,">");
	if (default_time) time_view=gsl_vector_subvector(default_time,j,1);
	time=sbtab_column_to_gsl_vector(DataTable_j,"!Time");
	assert(time);
	for (i=0;i<L;i++){
	  data_sub=gsl_matrix_submatrix(Y_DATA,i,0,1,Y_DATA->size2);
	  stdv_sub=gsl_matrix_submatrix(Y_STDV,i,0,1,Y_STDV->size2);
	  input_row=gsl_matrix_row(input_block,i);
	  time_view=gsl_vector_subvector(time,i,1);
	  map_index(IdxMap,j,i);
	  D=data_block_alloc(E_j,j,i,
			     &(data_sub.matrix),
			     &(stdv_sub.matrix),
			     &(input_row.vector),
			     &(time_view.vector),
			     lflag[j]);
	  g_ptr_array_add(Data,D);
	}
	break;
      default:
	fprintf(stderr,"[%s] error: unknown experiment type\n",__func__);
	abort();	
      }
    } else {
      fprintf(stderr,"Either (Level 1) Output or DataTable %i of %i missing (%p).\n",j,nE,DataTable_j);
      fflush(stdout);
      abort();
    }
  }
  return Data;
}

void write_data(gpointer data, gpointer user_data){
  data_t *D=data;
  data_attr_buffer_t *buffer=user_data;
  h5block_t *h5data=buffer->h5data;
  h5block_t *h5stdv=buffer->h5stdv;
  norm_t *N=buffer->N;
  int major=D->MajorIndex;
  int minor=D->MinorIndex;
  int lflag=D->lflag;
  int index=flat_index(buffer->IdxMap,major,minor);
  printf("[%s] writing data set with index %i (MAJOR=%i, MINOR=%i)\n",
	 __func__,index,major,minor);
  char *H5_data_name=D->Name->str;
  char *H5_stdv_name=D->Name->str; //same
  printf("[%s] creating dataset «%s».\n",__func__,H5_data_name);
  fflush(stdout);
  h5data->size[0]=D->measurement->size1;
  h5data->size[1]=D->measurement->size2;
  h5stdv->size[0]=D->noise->size1;
  h5stdv->size[1]=D->noise->size2;
  
  H5LTmake_dataset_double(h5data->group_id,H5_data_name,2,h5data->size,D->measurement->data);
  H5LTmake_dataset_double(h5stdv->group_id,H5_stdv_name,2,h5stdv->size,D->noise->data);
  H5LTset_attribute_int(h5data->group_id,H5_data_name,"LikelihoodFlag",&lflag, 1);
  H5LTset_attribute_int(h5data->group_id,H5_data_name,"index",&index, 1);
  H5LTset_attribute_int(h5data->group_id,H5_data_name,"major",&major, 1);
  H5LTset_attribute_int(h5data->group_id,H5_data_name,"minor",&minor, 1);
  H5LTset_attribute_double(h5data->group_id,H5_data_name,"time",D->time->data, D->time->size);
  H5LTset_attribute_double(h5data->group_id,H5_data_name,"input",D->input->data, D->input->size);
  /* write normalisation attributes */
  int rk, rt;
  if (N && N->experiment->len>0){
    if (norm_index(N->experiment,index,&rk)){
      H5LTset_attribute_int
	(h5data->group_id,H5_data_name,
	 "NormaliseByExperiment",&rk,1);
    }
  }
  if (N && N->time->len>0){
    if (norm_index(N->time,index,&rt)){
      H5LTset_attribute_int
	(h5data->group_id,H5_data_name,
	 "NormaliseByTimePoint",&rt,1);
    }
  }
  if (N && N->output->len>0){
    H5LTset_attribute_int
      (h5data->group_id,H5_data_name,
       "NormaliseByOutput",
       (int*) N->output->data,N->output->len);
  }
}

