#include "data.h"

data_t* data_block_alloc(int major, int minor, gsl_matrix *data, gsl_matrix *stdv, gsl_vector *input, gsl_vector *time){
  assert(time->size == data->size1);
  assert(data->size1 == stdv->size1 && data->size2 == stdv->size2);
  data_t *D=malloc(sizeof(data_t));
  D->MajorIndex=major;
  D->MinorIndex=minor;
  D->measurement=data;
  D->noise=stdv;
  D->input=input;
  D->time=time;
  return D;
}

GPtrArray* /* Data blocks. Each can be simulated with exactly one input vector. */
create_data
(GPtrArray *DataTable, /* Data tables exactly as they appear in the spreadsheet*/
 int *lflag, /* a flag that is true if the measurement contributes to Likelihood*/
 experiment_type *ExperimentType, /* Type: _Dose Response_ or _Time Series_ */
 GPtrArray *input_override, /* default input for each experiment */
 gsl_vector *default_time, /* default measurement time if not otherwise specified*/
 GArray **Normalisation, /* experiment normalisation structure, as three integer valued vectors */
 sbtab_t *Input,
 sbtab_t *Output)
{
  int j;
  int N,M,L;
  gsl_matrix Y_DATA, Y_STDV;
  gsl_vector *u_j, *time;
  int nE=DataTable->len;
  sbtab_t *DataTable_j;
  GPtrArray *Data=g_ptr_array_sized_new(nE);
  data_t *D;

  GPtrArray *uid=sbtab_find_column(Input,"!ID",NULL);
  GPtrArray *oid=sbtab_find_column(Output,"!ID",NULL);
  GPtrArray *ErrorName=sbtab_find_column(Output,"!ErrorName",NULL);
  gsl_matrix *input_block;
  gsl_vector_view input_row, time_row;
  gsl_matrix_view data_sub, stdv_sub;
  int nU=uid->len;
  
  for (j=0;j<nE;j++) { // j is the major experiment index
    DataTable_j=g_ptr_array_index(DataTable,j);
    u_j=g_ptr_array_index(input_override,j);
    if (DataTable_j){
      switch(ExperimentType[j]){
      case time_series:
	//	Y=get_data_matrix(DataTable_j,Output,lflag[j]);
	Y_DATA=sbtab_columns_to_gsl_matrix(DataTable_j, oid, ">", 1.0);
	Y_STDV=sbtab_columns_to_gsl_matrix(DataTable_j, ErrorName, NULL, INFINITY);
	time=sbtab_column_to_gsl_vector(DataTable_j,"!Time");
	D=data_block_alloc(j,0,Y_DATA,Y_STDV,u_j,time);
	g_ptr_array_append(Data,D);
	break;
      case dose_response:
	Y=get_data_matrix(DataTable_j,Output,lflag[j]);
	L=table_length(DataTable_j);
	input_block=gsl_matrix_alloc(L,nU);
	for (i=0;i<L;i++) {
	  input_row=gsl_matrix_row(input_block,i);
	  gsl_vector_memcpy(&(input_row.vector),u_j);
	}
	sbtab_update_gsl_matrix(input_block,DataTable_j,uid,">");
	if (default_time) time_view=gsl_vector_subvector(default_time,j,1);
	time=sbtab_column_to_gsl_vector(DataTable_j,"!Time");
	assert(time);
	N=Y[DATA]->size1;
	M=Y[DATA]->size2;
	assert(N>0 && M>0);
	for (i=0;i<N;i++){
	  data_sub=gsl_matrix_submatrix(Y[DATA],i,0,1,M);
	  stdv_sub=gsl_matrix_submatrix(Y[STDV],i,0,1,M);
	  input_row=gsl_matrix_row(input_block,i);
	  time_view=gsl_vector_subvector(time,i,1);
	  D=data_block_alloc(j,i,
			     &(data_sub.matrix),
			     &(sd_data_sub.matrix),
			     ,&(input_row.vector),
			     &(time_view.vector));
	  g_ptr_array_append(Data,D);
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

normalisation_t norm_alloc(int num_sim, int num_out){
  int n=size>0?size:10;
  normalisation_t *N=malloc(sizeof(normalisation_t));
  N->output=g_array_sized_new(FALSE,FALSE,sizeof(int),num_out);
  N->time=g_array_sized_new(FALSE,FALSE,sizeof(int),num_sim);
  N->experiment=g_array_sized_new(FALSE,FALSE,sizeof(int),num_sim);
  return N;
}


typedef struct {
  regex_t *re;
  sbtab_t *table;
  GPtrArray *sub_table;
  GArray **index;
  GArray *major;
  GArray *minor;
} user_data_t;

/* g_ptr_array_foreach function*/
void find_reference(gpointer data, gpointer user_data){
  gchar *field=data;
  user_data_t *buffer=user_data;

  gchar *Reference;
  int I=-1;
  int i=-1;
  sbtab_t *st;

  int re_n=(sub_table==NULL?1:2);
  
  if (field){
    GPtrArray *m=re_match(buffer->re,field,re_n);
    if (m->len>=2){
      Reference=g_ptr_array_index(m,1);
      I=sbtab_find_row(buffer->table->row,Reference->str)
    }
    if (I>=0 && sub_table && m->len>=3){
      st=g_ptr_array_index(sub_table,I);
      Reference=g_ptr_array_index(m,2);
      i=sbtab_find_row(st->row,Reference);
    }
    g_ptr_array_free(m);
  }
  g_array_add(buffer->major,I);
  g_array_add(buffer->minor,i);
}

normalisation_t* normalisation(sbtab_t *ExperimentTable, experiment_type *ExperimentType, sbtab_t *OutputTable, GPtrArray *DataTable, GArray **SimUnitIdx, GArray *MajorIdx, GArray *MinorIdx){
  regex_t ID_TimePoint, ID;
  GPtrArray *RelativeTo;
  int status;

  status=egrep(&ID_TimePoint,"([[:alpha:]][[:alnum:]_]*)\\[?([^]]+)?\\]?");
  status=egrep(&ID,"([[:alpha:]][[:alnum:]_]*)");
  
  gchar *field;
  int nE=table_length(ExperimentTable);
  guint K=MajorIdx->len;
  int nO=table_length(OutputTable);

  normalisation_t *N=norm_alloc(K,nO);
  user_data_t buffer;  
  assert(SimUnitIdx);
  assert(MajorIdx);
  assert(MinorIdx);

  RelativeTo=sbtab_find_column(OutputTable,"!RelativeTo !NormalisedBy !NormalizedBy",NULL);

  if (RelativeTo){
    buffer.sub_i=1;
    buffer.re=&ID;
    buffer.table=OutputTable;
    buffer.sub_table=NULL;
    buffer.major=N->output;
    g_ptr_array_foreach(RelativeTo,find_reference, &buffer);
  }
  /* experiment normalisation and time normalisation */
  RelativeTo=sbtab_find_column(ExperimentTable,"!RelativeTo !NormalisedBy !NormalizedBy",NULL);
  int I,i;
  int ri,rt;
  int rk,rI;
  if (RelativeTo){
    buffer.re=&ID_TimePoint;
    buffer.table=ExperimentTable;
    buffer.sub_table=DataTable;
    buffer.major=g_array_sized_new(FALSE,FALSE,sizeof(int),nE);
    buffer.minor=g_array_sized_new(FALSE,FALSE,sizeof(int),nE);
    g_ptr_array_foreach(RelativeTo,find_reference, &buffer);
    for (k=0;k<K;k++){
      I=g_array_index(MajorIdx,int,k);
      i=g_array_index(MinorIdx,int,k);
      ri=i;
      rt=-1;
      printf("[%s] Experiment %i.%i normalised by",__func__,I,i);
      fflush(stdout);
      switch (ExperimentType[I]){
      case dose_response:
	ri=g_array_index(buffer->minor,I);
	rt=0; 
	break;
      case time_series:
	ri=0;
	rt=minor[0];	  
	break;
      default:
	fprintf(stderr,"[%s] unknown experiment type: %i.\n",__func__,ExperimentType[I]);
	abort();
      }
      
      rk=-1;
      rI=g_array_index(buffer->major,int,I);
      assert(rI<nE);
      if (rI>=0 && SimUnitIdx[rI]!=NULL && ri<SimUnitIdx[rI]->len && ri>=0)
	rk=g_array_index(SimUnitIdx[rI],int,ri);
      g_array_append_val(N->experiment,rk);
      g_array_append_val(N->time,rt);      
    }
    fflush(stdout);
  }
  regfree(&ID_TimePoint);
  regfree(&ID);
  return N;
}
