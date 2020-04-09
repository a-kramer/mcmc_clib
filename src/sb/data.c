#include "data.h"

data_t* data_block_alloc(int major, int minor, gsl_matrix *data, gsl_matrix *stdv){
  data_t *D=malloc(sizeof(data_t));
  D->MajorIndex=major;
  D->MinorIndex=minor;
  D->measurement=data;
  D->noise=stdv;
  return D;
}

GPtrArray* create_data
(GPtrArray *DataTable, 
 int *lflag,
 experiment_type *ExperimentType,
 GPtrArray *E_default_input, /* default input for each experiment*/
 gsl_vector *default_time, /* default measurement time if not otherwise specified*/
 GArray **Normalisation, /* experiment normalisation structure, as three integer valued vectors */
 sbtab_t *Input,
 sbtab_t *Output)
{
  int j;
  gsl_matrix *Y;
  int nE=DataTable->len;
  sbtab_t *DataTable_j;
  GPtrArray *Data=g_ptr_array_sized_new(nE);
  data_t *D;

  input_ID=sbtab_get_column(Input,"!ID");

  for (j=0;j<nE;j++) { // j is the major experiment index
    DataTable_j=g_ptr_array_index(DataTable,j);
    if (DataTable_j){
      switch(ExperimentType[j]){
      case unknown_type:
	fprintf(stderr,"[%s] error: unknown experiment type\n",__func__);
	abort();
	break;
      case time_series:
	Y=get_data_matrix(DataTable_j,Output,lflag[j]);
	time=sbtab_column_to_gsl_vector(DataTable_j,"!Time");
	D=data_block_alloc(j,0,Y[DATA],Y[STDV]);
	break;
      case dose_response:
	Y=get_data_matrix(DataTable_j,Output,lflag[j]);
	input_block=get_input_matrix(DataTable_j,input_ID,E_default_input[j]);
	if (default_time) time_view=gsl_vector_subvector(default_time,j,1);
	time=sbtab_column_to_gsl_vector(DataTable_j,"!Time");
	N=Y[DATA]->size1;
	M=Y[DATA]->size2;
	assert(N>0 && M>0);
	for (i=0;i<N;i++){ // i is the minor "simulation unit" index
			   // within a dose response experiment
	  data_sub=gsl_matrix_submatrix(Y_dY[DATA],i,0,1,M);
	  sd_data_sub=gsl_matrix_submatrix(Y_dY[STDV],i,0,1,M);
	  input_row=gsl_matrix_row(input_block,i);
	  if (time) time_view=gsl_vector_subvector(time,i,1);
	  write_data_to_hdf5(file_id,
			     &(data_sub.matrix),
			     &(sd_data_sub.matrix),
			     &(time_view.vector),
			     &(input_row.vector),j,i,Normalisation,lflag[j]);	  
	}
	break;
      }	
      gsl_matrix_free(Y[DATA]);
      gsl_matrix_free(Y[STDV]);
      free(Y);
      if (time) gsl_vector_free(time);
    } else {
      fprintf(stderr,"Either (Level 1) Output or DataTable %i of %i missing (NULL).\n",j,nE);
      abort();
    }
  }



}
