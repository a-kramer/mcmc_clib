#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <glib.h>
#include <string.h>
#include <regex.h>
#ifndef _GNU_SOURCE
#include <libgen.h>
#endif
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "sbtab.h"
#include "../mcmc/model_parameters_smmala.h"
#include "../mcmc/ptype.h"
// "length" means: number of rows, for a table
#define get_table_length(Table) ((Table && Table->column) ? Table->column[0]->len:0)

#define DEFAULT_STR_LENGTH 128
// Experiment Types
#define UNKNOWN_TYPE 0
#define DOSE_RESPONSE 1
#define TIME_SERIES 2
// in an array, these are the indices for data and standrad deviation
#define DATA 0
#define STDV 1
// Normalisation is stored in a pointer array, using this order:
#define NORM_STRIDE 3
#define NORM_EXPERIMENT 0
#define NORM_OUTPUT 1
#define NORM_TIME 2

#define LOG_SCALE 0
#define LOG10_SCALE 1
#define LIN_SCALE 2


sbtab_t* parse_sb_tab(char *);
gsl_matrix** get_data_matrix(sbtab_t *DataTable, sbtab_t *L1_OUT, int IsNormalisingExperiment);
gsl_vector* get_time_vector(sbtab_t *DataTable);
gsl_matrix* get_input_matrix(sbtab_t *DataTable, GPtrArray *input_ID, gsl_vector *default_input);
herr_t write_data_to_hdf5(hid_t file_id, gsl_matrix *Y, gsl_matrix *dY, gsl_vector *time, gsl_vector *input, int major, int minor, GArray **N, int lflag);
herr_t write_prior_to_hdf5(hid_t file_id, gsl_vector *mu, void *S, int type);
int process_data_tables(hid_t file_id,  GPtrArray *sbtab,  GHashTable *sbtab_hash);
gchar* dup_match(regmatch_t *match, char *source);
gsl_vector* sbtab_column_to_gsl_vector(sbtab_t *table, gchar *column_name);
herr_t process_prior(hid_t file_id, GPtrArray *sbtab, GHashTable *sbtab_hash);

int needs_normalisation(GArray **N, int i){
  int by_experiment=N[NORM_EXPERIMENT] && g_array_index(N[NORM_EXPERIMENT],int,i)>=0;
  int by_timepoint=N[NORM_TIME] && g_array_index(N[NORM_TIME],int,i)>=0;
  int by_output=N[NORM_OUTPUT] && N[NORM_OUTPUT]->len>0;
  return (by_experiment || by_timepoint || by_output);
}

int main(int argc, char*argv[]){
  int i;
  int status=EXIT_SUCCESS;
  regex_t *sb_tab_file, *h5_file;
  sbtab_t *table;
  GPtrArray *sbtab;
  GHashTable *sbtab_hash;
  gchar *H5FileName=NULL;
  
  sb_tab_file=malloc(sizeof(regex_t));
  h5_file=malloc(sizeof(regex_t));
  regcomp(sb_tab_file,".*[.][tc]sv",REG_EXTENDED);
  regcomp(h5_file,".*[.]h5",REG_EXTENDED);
  sbtab=g_ptr_array_new_full(3,sbtab_free);
  sbtab_hash=g_hash_table_new(g_str_hash, g_str_equal);
  for (i=1;i<argc;i++){
    if (regexec(sb_tab_file,argv[i],0,NULL,0)==0) {
      printf("[main] parsing sbtab file %s.\n",argv[i]);
      table=parse_sb_tab(argv[i]);
      g_ptr_array_add(sbtab,table);
      g_hash_table_insert(sbtab_hash,table->TableName,table);
    } else if (regexec(h5_file,argv[i],0,NULL,0)==0) {
      H5FileName=g_strdup(argv[i]);
    } else printf("unknown option «%s».\n",argv[i]);
  }
  if (H5FileName==NULL){
      H5FileName=g_strdup("ExperimentalData.h5");
  }
  guint n_tab=sbtab->len;
  printf("[main] %i tables read: ",n_tab);
  for (i=0;i<n_tab;i++){
    table=g_ptr_array_index(sbtab,i);
    printf(" %s ",table->TableName);
  }
  printf("\n");
  hid_t file_id;
  file_id = H5Fcreate(H5FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  status &= process_data_tables(file_id, sbtab, sbtab_hash);
  status &= process_prior(file_id, sbtab, sbtab_hash);
  status &= H5Fclose(file_id);
  //process_input_table(H5FileName, sbtab, sbtab_hash);
  //process_prior(H5FileName, sbtab, sbtab_hash);
  return status;
}

gsl_vector* get_default_input(sbtab_t *Input){
  int j;
  double val;
  guint nU=get_table_length(Input);
  gsl_vector *default_input;
  default_input=gsl_vector_alloc(nU);
  gsl_vector_set_zero(default_input);
  GPtrArray *u_col;
  u_col=sbtab_get_column(Input,"!DefaultValue");
  printf("[get_default_input] ");
  if (u_col!=NULL){
    for (j=0;j<nU;j++){
      val=strtod(g_ptr_array_index(u_col,j),NULL);
      gsl_vector_set(default_input,j,val);
      printf(" %g ",val);
    }
    printf(";\n");
  } else {
    printf("none defined, defaulting to zero for all inputs.\n");
  }
  return default_input;
}

sbtab_t* find_input_table(GHashTable *sbtab_hash){
  sbtab_t *Input;
  Input=g_hash_table_lookup(sbtab_hash,"Input");
  assert(Input!=NULL);
  return Input;
}

sbtab_t* find_output_table(GHashTable *sbtab_hash){
  sbtab_t *L1_OUT;
  L1_OUT=g_hash_table_lookup(sbtab_hash,"L1_OUT");
  if (L1_OUT==NULL){
    L1_OUT=g_hash_table_lookup(sbtab_hash,"Output");
  }
  assert(L1_OUT!=NULL);
  return L1_OUT;
}

sbtab_t* find_experiment_table(GHashTable *sbtab_hash){
  sbtab_t *E_table;
  E_table=g_hash_table_lookup(sbtab_hash,"Experiments");
  if (E_table==NULL){
    E_table=g_hash_table_lookup(sbtab_hash,"TableOfExperiments");
  }
  assert(E_table!=NULL);
  return E_table;
}

gsl_vector *get_default_time(sbtab_t *E_table){
  gsl_vector *default_time;
  GPtrArray *default_time_column;
  default_time_column=sbtab_get_column(E_table,"!Time");
  int nE=get_table_length(E_table);
  int j;
  gchar *s;
  double t;
  default_time=gsl_vector_alloc(nE);
  printf("[get_default_time] ");
  if (default_time_column!=NULL){
    for (j=0;j<nE;j++) {
      s=g_ptr_array_index(default_time_column,j);
      t=strtod(s,NULL);
      gsl_vector_set(default_time,j,t);
      printf(" %g ",t);
    }
    printf(";\n");
  } else {    
    gsl_vector_set_zero(default_time);
    printf("not defined, setting to zero.\n");
  }
  return default_time;
}

int* get_experiment_type(GPtrArray *E_type){
  int j;
  int *ExperimentType;
  gchar *Type;
  regex_t DoseResponse;
  regcomp(&DoseResponse,"Dose[[:blank:]]*Response", REG_EXTENDED | REG_ICASE);
  regex_t TimeSeries;
  regcomp(&TimeSeries,"Time[[:blank:]]*Series", REG_EXTENDED | REG_ICASE);
  int nE=E_type->len;
  printf("[get_experiment_type] %i Experiments:",nE);
  ExperimentType=malloc(nE*sizeof(int));  
  for (j=0;j<nE;j++){
    Type=NULL;
    if (E_type!=NULL){
      Type=g_ptr_array_index(E_type,j);
      if (Type!=NULL){
	printf( "«%s»",Type);
	if (regexec(&DoseResponse, Type, 0, NULL, 0)==0){
	  ExperimentType[j]=DOSE_RESPONSE;
	} else if (regexec(&TimeSeries, Type, 0, NULL, 0)==0){
	  ExperimentType[j]=TIME_SERIES;  
	} else {
	  //fprintf(stderr,"[get_experiment_type] This type is not known: «%s».\n",Type);
	  ExperimentType[j]=UNKNOWN_TYPE;
	}
      } else {
	ExperimentType[j]=UNKNOWN_TYPE;
	fprintf(stderr,"[get_experiment_type] Experiment Type is set, but is not one of: «Time Series», «Dose Response»).\n");	
      }      
    }else{
      // default case, when !Type has not been specified at all, for now
      ExperimentType[j]=TIME_SERIES;
      printf("[get_experiment_type] Experiment Type is not defined, assuming time series data. Use !Type column in Table of Experiments (Possible Values: «Time Series», «Dose Response»).\n");
    }
    printf("(%i) ",ExperimentType[j]);
  }
  printf(";\n");
  regfree(&DoseResponse);
  regfree(&TimeSeries);
  fflush(stdout);
  return ExperimentType;
}

gsl_vector** get_experiment_specific_inputs(sbtab_t *ExperimentTable, sbtab_t *Input, gsl_vector *default_input){
  int j,k;
  gchar *s, *input_id;
  GPtrArray *u_col;
  double val;
  GPtrArray *input_ID;
  //gchar *experiment_id;
  //gchar *experiment_name;
  gsl_vector **E_default_input;
  int nE=get_table_length(ExperimentTable);
  E_default_input=malloc(sizeof(gsl_vector*)*nE);
  assert(default_input!=NULL);
  int nU=default_input->size;
  for (j=0;j<nE;j++) { // j is the major experiment index
    E_default_input[j]=gsl_vector_alloc(nU);
    gsl_vector_memcpy(E_default_input[j],default_input);
    
    input_ID=sbtab_get_column(Input,"!ID");
    assert(input_ID!=NULL);
    
    for (k=0;k<nU;k++){
      input_id=g_ptr_array_index(input_ID,k);
      u_col=sbtab_get_column(ExperimentTable,input_id);
      if (u_col!=NULL){
	s=g_ptr_array_index(u_col,j);
	val=strtod(s,NULL);
	gsl_vector_set(E_default_input[j],k,val);
      }
    }
  }
  return E_default_input;
}

GArray** get_normalisation(sbtab_t *ExperimentTable, int *ExperimentType, sbtab_t *OutputTable, sbtab_t **DataTable, GArray **SimUnitIdx, GArray *MajorIdx, GArray *MinorIdx){
  regex_t ID_TimePoint, ID;
  GPtrArray *RelativeTo;
  int status;
  char errbuf[DEFAULT_STR_LENGTH];
  status=regcomp(&ID_TimePoint,"([[:alpha:]][[:alnum:]_]*)\\[?([^]]+)?\\]?", REG_EXTENDED);
  if (status!=0){
    regerror(status,&ID_TimePoint,errbuf,DEFAULT_STR_LENGTH);
    perror(errbuf);
  }
  status=regcomp(&ID,"([[:alpha:]][[:alnum:]_]*)", REG_EXTENDED);
  if (status!=0){
    regerror(status,&ID,errbuf,DEFAULT_STR_LENGTH);
    perror(errbuf);
  }
  guint *major, *minor;
  int i,k,I,rI,ri,rk,rt;
  int *rO;
  int j,l;
  gchar *field;
  gchar *RefExp;
  gchar *RefLine;
  gchar *RefOut;
  regmatch_t match[4];  
  int nE=get_table_length(ExperimentTable);
  guint K=MajorIdx->len;
  GArray **N;
  int nO=get_table_length(OutputTable);

  N=malloc(sizeof(GArray*)*NORM_STRIDE);

  assert(SimUnitIdx!=NULL);
  assert(MajorIdx!=NULL);
  assert(MinorIdx!=NULL);
  // output normalisation
  printf("[get_normalisation] %i outputs, normalised by: ",nO);
  RelativeTo=sbtab_get_column(OutputTable,"!RelativeTo");
  if (RelativeTo==NULL) RelativeTo=sbtab_get_column(OutputTable,"!NormalisedBy");
  fflush(stdout);
  if (RelativeTo!=NULL){
    N[NORM_OUTPUT]=g_array_sized_new(FALSE,FALSE,sizeof(int),nO);    
    for (j=0;j<nO;j++){      
      field=g_ptr_array_index(RelativeTo,j);
      assert(field!=NULL);
      printf(" «%s» ",field); fflush(stdout);
      status=regexec(&ID, field, 4, match, REG_EXTENDED);
      //printf("[%i]",status); fflush(stdout);
      if (status==0){
	RefOut=dup_match(&match[1], field);
	//RefOut=field;
	assert(RefOut!=NULL);
	printf("(%s=",RefOut);fflush(stdout);
	rO=g_hash_table_lookup(OutputTable->row,RefOut);
	l=(rO==NULL?(-1):rO[0]);
	printf("%i)",l); fflush(stdout);
	g_array_append_val(N[NORM_OUTPUT],l);
	g_free(RefOut);
      } else {
	fprintf(stderr,"[get_normalisation] field «%s» is not an ID.\n",field);
	l=-1;
	g_array_append_val(N[NORM_OUTPUT],l);
      }
      field=NULL;
    }
  } else {
    N[NORM_OUTPUT]=NULL;
  }
  printf("\n");
  fflush(stdout);
  // experiment normalisation and time normalisation
  printf("[get_normalisation] %i simulation units.\n",K); fflush(stdout);
  RelativeTo=sbtab_get_column(ExperimentTable,"!RelativeTo");
  if (RelativeTo==NULL) RelativeTo=sbtab_get_column(ExperimentTable,"!NormalisedBy");

  if (RelativeTo){
    N[NORM_EXPERIMENT]=g_array_sized_new(FALSE,FALSE,sizeof(int),K);
    N[NORM_TIME]=g_array_sized_new(FALSE,FALSE,sizeof(int),K);
    for (k=0;k<K;k++){  // normalisation array has to be of the same length as the number of simulation units
      I=g_array_index(MajorIdx,int,k);
      i=g_array_index(MinorIdx,int,k);
      printf("[get_normalisation] Experiment %i.%i normalised by",I,i);
      fflush(stdout);
      field=g_ptr_array_index(RelativeTo,I);
      regexec(&ID_TimePoint, field, 4, match, REG_EXTENDED);
      RefExp=dup_match(&match[1], field);
      RefLine=dup_match(&match[2], field);
      if (RefExp!=NULL && strlen(RefExp)>0){
	printf(" «%s»",RefExp);
	if (RefLine!=NULL) printf(" at «%s».\n",RefLine);
	else printf(" (point by point).\n");
	major=g_hash_table_lookup(ExperimentTable->row,RefExp);
	assert(major!=NULL);
	rI=major[0];
      } else {
	printf(" ...none.\n");
	rI=-1;
      }
      minor=NULL;
      if (rI>0){
	if (RefLine){
	  if (DataTable && rI<nE){
	    minor=g_hash_table_lookup(DataTable[rI]->row,RefLine);
	  } else {
	    fprintf(stderr,"[get_normalisation] error: DataTable[%i] does not exist.\n",rI);
	    exit(-1);
	  }
	  // minor can refer to a TimePoint, a Dose, or not exist at all (because of typos).
	  
	  if (minor==NULL){
	    fprintf(stderr,"[get_normalisation] Experiment %i.%i has an invalid normalisation key (SBtab ID reference: «%s»).\n",I,i,field);
	    exit(-1);
	  } else if (ExperimentType[I]==DOSE_RESPONSE){
	    // then there's only one measurement time
	    ri=minor[0];
	    rt=0; 
	  } else if (ExperimentType[I]==TIME_SERIES){
	    // then there's only one simulation unit
	    ri=0;
	    rt=minor[0];	  
	  } else {
	    fprintf(stderr,"[get_normalisation] unknown experiment type: %i.\n",ExperimentType[I]);
	    ri=0;
	    rt=minor[0];	    
	  }
	} else {
	  // point by point normalisation
	  assert(DataTable[I]->column[0]->len == DataTable[rI]->column[0]->len);
	  assert(ExperimentType[I]==ExperimentType[rI]);
	  ri=i;
	  rt=-1;	
	}
	rk=g_array_index(SimUnitIdx[rI],int,ri);      
      }else {
	rk=-1;
      }
      g_array_append_val(N[NORM_EXPERIMENT],rk);
      g_array_append_val(N[NORM_TIME],rt);      
    }
    fflush(stdout);
  }else{
    N[NORM_EXPERIMENT]=NULL;
    N[NORM_TIME]=NULL;
  }
  regfree(&ID_TimePoint);
  regfree(&ID);
  return N;
}

sbtab_t* get_data_table(GHashTable *sbtab_hash, gchar *experiment_id, gchar *experiment_name){
  sbtab_t *DataTable;
  DataTable=g_hash_table_lookup(sbtab_hash,experiment_id);
  if (DataTable==NULL) {
    fprintf(stderr,"\t\tCould not find Measurement Table by ID:\t«TableName='%s'», \n\t\twill try searching by Name instead:\t«TableName='%s'».\n",experiment_id,experiment_name);
    if (experiment_name!=NULL) DataTable=g_hash_table_lookup(sbtab_hash,experiment_name);
  }
  assert(DataTable!=NULL);
  return DataTable;
}

/* This is a generic function to convert a column to a numeric vector
 * Mor specific variants of this exist for some Quantities such as
 * time-vectors and data-matrices. This one allows missing values.
 *
 * gsl_vector* sbtab_column_to_gsl_vector(gchar *column_name, sbtab_t *table)
 */
gsl_vector* sbtab_column_to_gsl_vector(sbtab_t *table,gchar *column_name){
  int j,n;
  GPtrArray *column;
  gsl_vector *v;
  gchar *s,*r;
  double value;
  column=sbtab_get_column(table,column_name);
  if (column!=NULL){
    n=column->len;
    v=gsl_vector_alloc(n);
    for (j=0;j<n;j++){
      s=g_ptr_array_index(column,j);
      r=s;
      value=strtod(s,&r);
      if (s==r) value=NAN;
      gsl_vector_set(v,j,value);
    }
  }else{
    v=NULL;
  }
  return v;
}

int  get_experiment_index_mapping(int nE, int *ExperimentType, sbtab_t **DataTable, GArray **SU_index, GArray *MajorExpIdx, GArray *MinorExpIdx){
  int k=0,j,i;
  int N; // number of dose response rows

  for (j=0;j<nE;j++) { // j is the major experiment index
    N=DataTable[j]->column[0]->len;
    if (DataTable!=NULL){
      switch(ExperimentType[j]){
      case UNKNOWN_TYPE:
	fprintf(stderr,"[get_experiment_index_mapping] error: unknown experiment type\n");
	exit(-1);
	break;
      case TIME_SERIES:
	i=0;
	g_array_append_val(MajorExpIdx,j);	
	g_array_append_val(MinorExpIdx,i);
	SU_index[j]=g_array_new(FALSE, TRUE, sizeof(int));
	g_array_append_val(SU_index[j],k);
	k++;
	break;
      case DOSE_RESPONSE:
	SU_index[j]=g_array_sized_new(FALSE, TRUE, sizeof(int), N);
	for (i=0;i<N;i++){ // i is the minor "simulation unit" index
			   // within a dose response experiment
	  g_array_append_val(MajorExpIdx,j);
	  g_array_append_val(MinorExpIdx,i);
	  g_array_append_val(SU_index[j],k);
	  k++;
	}
	break;
      }	
    }
  }
  return k;
}

int process_data_tables(hid_t file_id,  GPtrArray *sbtab,  GHashTable *sbtab_hash){
  sbtab_t *E_table, *Input;
  GPtrArray *LikelihoodFlag;  
  gchar *experiment_id;
  //gchar *input_id;
  gchar *experiment_name;
  //gsl_vector *input;
  gsl_matrix **Y_dY=NULL; // to make one return pointer possible: Y_dY[0] is data, Y_dY[1] is standard deviation
  sbtab_t **DataTable;
  sbtab_t *L1_OUT;
  GPtrArray *ID, *E_type, *input_ID, *E_Name; 
  int status=EXIT_SUCCESS;
  int i,j; // loop counters
  regex_t TrueRE, FalseRE;
  regcomp(&TrueRE,"true|1",REG_EXTENDED|REG_ICASE|REG_NOSUB);
  regcomp(&FalseRE,"false|0",REG_EXTENDED|REG_ICASE|REG_NOSUB);
  // get table pointers
  Input=find_input_table(sbtab_hash);
  L1_OUT=find_output_table(sbtab_hash);
  E_table=find_experiment_table(sbtab_hash);
  // hdf5 file
  hid_t data_group_id, sd_group_id;  
  
  data_group_id = H5Gcreate2(file_id, "/data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  sd_group_id = H5Gcreate2(file_id, "/stdv", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert(data_group_id>0 && sd_group_id>0);
  guint nE=get_table_length(E_table);
  printf("[process_data_tables] %i Experiments.\n",nE);
  ID=sbtab_get_column(E_table,"!ID");
  assert(ID!=NULL);
  if (ID!=E_table->column[0]) printf("[process_data_tables] warning: !ID is not first column, contrary to SBtab specification.\n");
  E_Name=sbtab_get_column(E_table,"!Name");
  E_type=sbtab_get_column(E_table,"!Type");
  LikelihoodFlag=sbtab_get_column(E_table,"!Likelihood");
  
  gsl_vector *time=NULL, *default_time=NULL;
  default_time=get_default_time(E_table);
  int *ExperimentType;  
  ExperimentType=get_experiment_type(E_type);
  
  printf("[process_data_tables] found list of experiments: %s with %i entries\n",E_table->TableName,nE);
  gsl_vector *default_input;
  default_input=get_default_input(Input);
  // each experiment can override the default input:
  gsl_vector **E_default_input;
  E_default_input=get_experiment_specific_inputs(E_table, Input, default_input);

  //int nU=default_input->size;
  //input=gsl_vector_alloc(nU);
  gsl_matrix *input_block;
  gsl_matrix_view data_sub, sd_data_sub;
  gsl_vector_view input_row;
  gsl_vector_view time_view;
  int N,M;
  char *flag;
  int lflag[nE];
  DataTable=malloc(sizeof(sbtab_t*)*nE);
  for (j=0;j<nE;j++) { // j is the major experiment index
    if (ID!=NULL) experiment_id=(gchar*) g_ptr_array_index(ID,j);
    else experiment_id=(gchar*) g_ptr_array_index(E_table->column[0],j);
      
    if (E_Name!=NULL) {
      experiment_name=(gchar*) g_ptr_array_index(E_Name,j);
    } else {
      experiment_name=NULL;
    }
    printf("[process_data_tables] processing %s\n",experiment_id);
    DataTable[j]=get_data_table(sbtab_hash, experiment_id, experiment_name);
    if (LikelihoodFlag){
      flag=g_ptr_array_index(LikelihoodFlag,j);
      lflag[j]=(regexec(&FalseRE,flag,0,NULL,0)!=0);
      if (!lflag[j]){
	printf("[process_data_tables] Experiment[%i] %s «%s» will be used for normalisation purposes only; it will not explicitly appear in the likelihood.\n",j,experiment_id,experiment_name);}
    } else {
      lflag[j]=1; // default case; !Likelihood field can deactivate evaluation
    }
  }
  regfree(&FalseRE);
  regfree(&TrueRE);
  printf("[process_data_tables] Likelihood contribution flag for all(%i) experiments:",nE);
  for (j=0;j<nE;j++) printf(" %i ",lflag[j]);
  printf(".\n");
  // These two arrays map the overall index k of "simulation units" to major and minor
  GArray *MajorExpIdx;
  GArray *MinorExpIdx;
  MajorExpIdx=g_array_sized_new(FALSE, TRUE, sizeof(int), nE);
  MinorExpIdx=g_array_new(FALSE, TRUE, sizeof(int));

  // this entity stores the simulation unit index k, given major and minor version
  GArray **SU_index;
  SU_index=malloc(sizeof(GArray*)*nE);  
  int NumSimUnits=get_experiment_index_mapping(nE,ExperimentType,DataTable,SU_index,MajorExpIdx,MinorExpIdx);
  printf("[process_data_tables] the model will have to be simulated %i times.\n",NumSimUnits);
  sbtab_t *OutputTable;
  OutputTable = g_hash_table_lookup(sbtab_hash,"Output");
  GArray **Normalisation;
  Normalisation=get_normalisation(E_table, ExperimentType, OutputTable, DataTable, SU_index, MajorExpIdx, MinorExpIdx);
  for (j=0;j<nE;j++) { // j is the major experiment index
    if (DataTable[j]!=NULL){
      switch(ExperimentType[j]){
      case UNKNOWN_TYPE:
	fprintf(stderr,"[process_data_table] error: unknown experiment type\n");
	exit(-1);
	break;
      case TIME_SERIES:
	Y_dY=get_data_matrix(DataTable[j],L1_OUT,lflag[j]);
	time=get_time_vector(DataTable[j]);
	write_data_to_hdf5(file_id,
			   Y_dY[DATA],
			   Y_dY[STDV],
			   time,
			   E_default_input[j],j,0,Normalisation,lflag[j]);
	break;
      case DOSE_RESPONSE:
	Y_dY=get_data_matrix(DataTable[j],L1_OUT,lflag[j]);
	input_ID=sbtab_get_column(Input,"!ID");
	input_block=get_input_matrix(DataTable[j],input_ID,E_default_input[j]);
	if (default_time) time_view=gsl_vector_subvector(default_time,j,1);
	time=get_time_vector(DataTable[j]);
	N=Y_dY[DATA]->size1;
	M=Y_dY[DATA]->size2;
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
      gsl_matrix_free(Y_dY[DATA]);
      gsl_matrix_free(Y_dY[STDV]);
      if (time) gsl_vector_free(time);
    } else {
      fprintf(stderr,"Either L1 Output or DataTable %s missing (NULL).\n",experiment_id);
      status&=EXIT_FAILURE;
    }
  }
  //printf("result of reading matrices (Y,dY):\n");
  //for (j=0;j<nE;j++) gsl_matrix_fprintf(stdout,Y[j],"%g,");
    
  status &= H5Gclose(data_group_id);
  status &= H5Gclose(sd_group_id);
  return status;
}

int make_regular_expression(regex_t *RE, char *pattern){
  int EC=0;
  char error_buffer[128];
  regcomp(RE,pattern,REG_EXTENDED|REG_ICASE);
  if (EC!=0) {
    regerror(EC,RE,error_buffer,128);
    perror(error_buffer);
  }  
  return EC;
}

gsl_vector_int* get_parameter_scale(GPtrArray *Scale){
  guint n=0;
  guint i;
  gsl_vector_int *scale_type=NULL;
  regex_t LogType;
  regex_t Log10Type;
  regex_t LinType;
  gchar *Type;
  regmatch_t match[5];
  int EC=0;
  assert(Scale);
  n=Scale->len;
  // int regexec(const regex_t *restrict preg, const char *restrict string,
  //             size_t nmatch, regmatch_t pmatch[restrict], int eflags);
  EC&=make_regular_expression(&LogType,"^(natural|base-e)?[[:blank:]]*(log(arithm)?)$|^ln$");
  EC&=make_regular_expression(&Log10Type,"^(decadic|base-10)[[:blank:]]*(log(arithm)?)$|^log10$");
  EC&=make_regular_expression(&LinType,"^lin(ear)$");
  assert(EC==0);
 
  scale_type=gsl_vector_int_alloc((size_t) n);
  gsl_vector_int_set_all(scale_type,-1);
  for (i=0;i<n;i++){
    Type=g_ptr_array_index(Scale,i);
    if (regexec(&LogType, Type, 3, match, 0)==0){
      //printf("[get_parameter_scale] (%i) «%s» matches RE for 'natural logarithm'\n",i,Type);
      gsl_vector_int_set(scale_type,i,LOG_SCALE);
    } else if (regexec(&Log10Type, Type, 3, match, 0)==0){
      //printf("[get_parameter_scale] (%i) «%s» matches RE for 'base-10 logarithm'\n",i,Type);
      gsl_vector_int_set(scale_type,i,LOG10_SCALE);
    } else if (regexec(&LinType, Type, 3, match, 0)==0){
      //printf("[get_parameter_scale] (%i) «%s» matches RE for 'linear scale'\n",i,Type);
      gsl_vector_int_set(scale_type,i,LIN_SCALE);
    } else {
      printf("[get_parameter_scale] This «!Scale[%i]» is unknown: «%s»\n",(int) i, Type);
      exit(-1);
    }
  }
  printf("[get_parameter_scale] ScaleType legend: natural logarithm = %i, base-10 logarithm = %i, linear = %i\n",LOG_SCALE,LOG10_SCALE,LIN_SCALE);
  for (i=0; i<n; i++) printf("%1i",gsl_vector_int_get(scale_type,i));
  printf("\n");
  //  gsl_vector_int_fprintf(stdout,scale_type,"%i");
  return scale_type;
}

/* the mcmc sampling software operates in natural logarithm scale
 * so everything but log-scale needs to be conerted to it
 */
int adjust_scale(gsl_vector *mu, gsl_vector *stdv, gsl_vector_int *ScaleType){
  assert(ScaleType);
  assert(mu);  
  int i,n=ScaleType->size;
  int t_i;
  double LOG10=gsl_sf_log(10);
  double mu_i, stdv_i, median_i;
  double es2; // exp(sigma^2) of logspace sigma
  for (i=0;i<n;i++){
    t_i=gsl_vector_int_get(ScaleType,i);
    mu_i=gsl_vector_get(mu,i);
    stdv_i=gsl_vector_get(stdv,i);
    switch (t_i){
    case LIN_SCALE:      
      median_i=mu_i;
      //var_i=gsl_pow_2(stdv_i);
      // find mu and sigma of logspace
      mu_i=gsl_sf_log(median_i);
      es2=0.5+sqrt(0.25+gsl_pow_2(stdv_i/median_i));
      assert(es2>0);
      stdv_i=sqrt(gsl_sf_log(es2));
      // update the values of mu and stdv, now in logspace
      gsl_vector_set(stdv,i,stdv_i);
      gsl_vector_set(mu,i,mu_i);
      break;    
    case LOG10_SCALE:
      // here we just have to do a very simple shift by a factor of log(10).
      mu_i*=LOG10;
      stdv_i*=LOG10;
      gsl_vector_set(mu,i,mu_i);
      gsl_vector_set(stdv,i,stdv_i);
      break;
    /* case LOG_SCALE: */
    /*   printf("[adjust_scale] (%i) prior value (parameter DefaultValues) seems to be in log-scale (natural logarithm); no scale conversion needed.\n",i); */
    }
  }
  return EXIT_SUCCESS;
}

herr_t process_prior(hid_t file_id, GPtrArray *sbtab, GHashTable *sbtab_hash){
  int i,n;
  printf("[process_prior] Looking up «Parameter» Table.\n");
  sbtab_t *Parameters;
  Parameters=g_hash_table_lookup(sbtab_hash,"Parameter");
  assert(Parameters!=NULL);
  GPtrArray *P_ID=sbtab_get_column(Parameters,"!ID");
  if (P_ID==NULL) P_ID=sbtab_get_column(Parameters,"!Name");
  
  assert(P_ID!=NULL);
  printf("[process_prior] Found an ID or Name column.\n");
  
  sbtab_t *Covariance, *Precision;
  Covariance=sbtab_find(sbtab_hash,"ParameterCovariance Covariance Sigma SIGMA cov");
  Precision=sbtab_find(sbtab_hash,"ParameterPrecision PriorPrecision Precision precision InvCovariance InverseCovariance inv(Sigma) inv(SIGMA)");
  printf("[process_prior] Checking Prior Type.\n");
  sbtab_t *cov[2]={Covariance, Precision};
  int any_cov_given=0;
  for (i=0;i<2;i++) any_cov_given+=(cov[i]!=NULL);
  printf("[process_prior] Any type of Covariance or Precision matrix given: %i.\n",any_cov_given); fflush(stdout);
  gchar *s;
  
  gsl_vector *mu=NULL;  
  gchar **ValueName;

  // get the scale (logarithmic or not)
  GPtrArray *Scale;
  Scale=sbtab_get_column(Parameters,"!Scale");
  gsl_vector_int *ScaleType; // 0 linear; 1 natural logarithm; 10 decadic logarithm;
  if (Scale){
    fprintf(stdout,"[process_prior] «!Scale» found in Parameter Table.\n");
    ScaleType=get_parameter_scale(Scale);
  }else{
    fprintf(stderr,"[process_prior] «!Scale» not found in Parameter Table, defaulting to LOG10 Scale.\n");
    ScaleType=gsl_vector_int_alloc((int) P_ID->len);
    gsl_vector_int_set_all(ScaleType,1);
  }
  
  mu=sbtab_column_to_gsl_vector(Parameters,"!DefaultValue");
  if (!mu){ // the above name is not a column, so try a couple of other names
    ValueName=g_strsplit("!Value !Mean !Median !Mode"," ",-1);  
    n=(int) g_strv_length(ValueName);
    i=-1;
    while (!mu && i<n){
      mu=sbtab_column_to_gsl_vector(Parameters,ValueName[++i]);
    }
    assert(mu);
    printf("[process_prior] Using column «%s» for prior μ (mu).\n",ValueName[i]);
    g_strfreev(ValueName);
  }
  int D=mu->size;

  void *S; // S is either Precision, Covariance, or sigma;
  int type=PRIOR_IS_UNKNOWN;
  gsl_matrix *M=NULL;
  printf("[process_prior] prior has size %i.\n",D);

  if (any_cov_given){
    sbtab_t *C=NULL;
    type=PRIOR_IS_GAUSSIAN;
    ALSO(type,PRIOR_IS_MULTIVARIATE);

    M=gsl_matrix_alloc(D,D);
    gsl_matrix_set_identity(M);
    gsl_vector_view column;
    gsl_vector *c;
    if (Precision!=NULL){
      ALSO(type,PRIOR_PRECISION_GIVEN);
      C=Precision;
    } else if (Covariance!=NULL) {
      ALSO(type,PRIOR_SIGMA_GIVEN);
      C=Covariance;
    } else {
      printf("[process_prior] this prior specification {%i} is not handled yet.\n",any_cov_given);
      exit(-1);
    }
    // regardless of inversion issues, read all the values into the gsl_matrix
    assert(C);
    for(i=0;i<D;i++){
      column=gsl_matrix_column(M,i);
      s=g_ptr_array_index(P_ID,i);
      c=sbtab_column_to_gsl_vector(C,s);
      if (c){
    	gsl_vector_memcpy(&(column.vector),c);
	gsl_vector_free(c);
      }
    }
    S=M;
  } else {
    gsl_vector *sigma;
    sigma=gsl_vector_alloc(D);
    type=PRIOR_IS_GAUSSIAN;
    gsl_vector *stdv;
    stdv=sbtab_column_to_gsl_vector(Parameters,"!Std");
    if (stdv){
      printf("[process_prior] using !Std column.\n");
      adjust_scale(mu,stdv,ScaleType);
      gsl_vector_memcpy(sigma,stdv);
      gsl_vector_free(stdv);      
      S=sigma;
    } else {
      ALSO(type,PRIOR_IS_GENERALISED);
      gsl_vector *par_max, *par_min;
      gsl_vector *alpha;
      par_max=sbtab_column_to_gsl_vector(Parameters,"!Max");
      par_min=sbtab_column_to_gsl_vector(Parameters,"!Min");
      if (par_max && par_min){
	alpha=gsl_vector_alloc(D);
	gsl_vector_set_zero(alpha);
	printf("[process_prior] Creating generalised Gaussian prior from !Min and !Max fields.\n");
	// apha is the width of the [min,max] interval 
	gsl_vector_add(alpha,par_max);
	gsl_vector_sub(alpha,par_min);
	gsl_vector_scale(alpha,0.5);
	S=alpha;
      }else{
	printf("[process_prior] neither !Std, nor (!Max,!Min) pairs were specified.\n");
	exit(-1);
      }
    }     
  }
  herr_t status = write_prior_to_hdf5(file_id, mu, S, type);
  return status;
}

herr_t write_prior_to_hdf5(hid_t file_id, gsl_vector *mu, void *S, int type){
  herr_t status=0;
  hid_t prior_group_id;
  hsize_t mu_size=mu->size;
  gsl_matrix *Sigma;
  gsl_vector *sigma;
  gchar *Name=NULL;
  prior_group_id=H5Gcreate(file_id,"/prior",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status &= H5LTmake_dataset_double(prior_group_id,"mu",1,&mu_size,mu->data);
  printf("[write_prior_to_hdf5] prior type %i\n",type);
  if (PTYPE(type,PRIOR_IS_GAUSSIAN)){    
    if (PTYPE(type,PRIOR_IS_MULTIVARIATE)){
      printf("[write_prior_to_hdf5] prior is multivariate Gaussian.\n"); fflush(stdout);
      Sigma=(gsl_matrix*) S;      
      hsize_t Sigma_size[2];
      Sigma_size[0]=Sigma->size1;
      Sigma_size[1]=Sigma->size2;
      if (PTYPE(type,PRIOR_SIGMA_GIVEN))           Name=g_strdup("Sigma");
      else if (PTYPE(type,PRIOR_PRECISION_GIVEN))  Name=g_strdup("Precision");
      status &= H5LTmake_dataset_double(prior_group_id,Name,2,Sigma_size,Sigma->data);
    } else {
      printf("[write_prior_to_hdf5] prior is product of Gaussians.\n"); fflush(stdout);
      sigma=(gsl_vector*) S;
      if (PTYPE(type,PRIOR_IS_GENERALISED))        Name=g_strdup("alpha");
      else                                         Name=g_strdup("sigma");
      hsize_t sigma_size[1];
      sigma_size[0]=(hsize_t) sigma->size;
      status &= H5LTmake_dataset_double(prior_group_id,Name,1,sigma_size,sigma->data);
    }
    assert(prior_group_id>=0);
    printf("[write_prior_to_hdf5] prior_group_id = %li.\n",prior_group_id);
    //status &= H5Gopen(file_id,"/prior",H5P_DEFAULT);
  } else {
    printf("[write_prior_to_hdf5] Case not handled: %i\n",type);
    exit(-1);
  }
  if (Name) g_free(Name);
  status &= H5Gclose(prior_group_id);
  return status;
}

herr_t write_data_to_hdf5(hid_t file_id, gsl_matrix *Y, gsl_matrix *dY, gsl_vector *time, gsl_vector *input, int major, int minor, GArray **N, int lflag){
  herr_t status;
  hid_t data_group_id, sd_group_id, dataspace_id, dataset_id, sd_data_id;  
  hsize_t size[2];
  char H5_data_name[128];
  static int index=0;
  /* Create a new file using default properties. */
  status=EXIT_SUCCESS;
  printf("[write_data_to_hdf5] writing data set with index %i (MAJOR=%i, MINOR=%i)\n",index,major,minor);
  printf("[write_data_to_hdf5] looking up data group «H5_ROOT/data»");
  data_group_id=H5Gopen(file_id,"/data",H5P_DEFAULT); printf(" id=%li\n",data_group_id);
  printf("[write_data_to_hdf5] looking up standard deviation group «H5_ROOT/stdv»");
  sd_group_id=H5Gopen(file_id,"/stdv",H5P_DEFAULT); printf(" id=%li\n",sd_group_id);
  assert(data_group_id>0);
  assert(sd_group_id>0);
  // write data and standard deviation to file
  size[0]=Y->size1;
  size[1]=Y->size2;
  dataspace_id = H5Screate_simple(2, size, NULL);
  // Y
  sprintf(H5_data_name,"data_block_%i",index);
  printf("[write_data_to_hdf5] creating dataset «%s».\n",H5_data_name); fflush(stdout);
  dataset_id = H5Dcreate2(data_group_id, H5_data_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  assert(dataset_id>0);
  status &= H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Y->data);
  status &= H5LTset_attribute_int(data_group_id,H5_data_name,"LikelihoodFlag",&lflag, 1);
  status &= H5LTset_attribute_int(data_group_id,H5_data_name,"index",&index, 1);
  status &= H5LTset_attribute_int(data_group_id,H5_data_name,"major",&major, 1); // major experiment index, as presented in SBtab file
  status &= H5LTset_attribute_int(data_group_id,H5_data_name,"minor",&minor, 1); // minor experiment index, within a dose response experiment
  // time
  status &= H5LTset_attribute_double(data_group_id,H5_data_name,"time",time->data, time->size);
  // input
  status &= H5LTset_attribute_double(data_group_id,H5_data_name,"input",input->data, input->size);
  // Normalisation attributes
  int *rK,*rT,*rO;
  if (N!=NULL){
    if (N[NORM_EXPERIMENT]!=NULL && N[NORM_EXPERIMENT]->len>index){
      rK=&g_array_index(N[NORM_EXPERIMENT],int,index);
      //printf("[write_data] Simulation Unit %i is normalised by unit %i .. ",index,rK[0]);
      if (rK[0]>=0){
	//printf("creating");
	status &= H5LTset_attribute_int(data_group_id,H5_data_name,"NormaliseByExperiment",rK,1);
      } //else {
	//printf("omitting");
      //}
      //printf(" attribute.\n");
    } 
    if (N[NORM_TIME]!=NULL && N[NORM_TIME]->len>index){
      rT=&g_array_index(N[NORM_TIME],int,index);
      if (rT[0]>0){
	status &= H5LTset_attribute_int(data_group_id,H5_data_name,"NormaliseByTimePoint",rT,1);
      }
    }
    if (N[NORM_OUTPUT]!=NULL && N[NORM_OUTPUT]->len>0){
      rO=&g_array_index(N[NORM_OUTPUT],int,0);
      status &= H5LTset_attribute_int(data_group_id,H5_data_name,"NormaliseByOutput",rO,N[NORM_OUTPUT]->len);
    }    
  }
  status &= H5Dclose(dataset_id);
  // dY
  sprintf(H5_data_name,"stdv_block_%i",index);
  sd_data_id = H5Dcreate2(sd_group_id, H5_data_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status &= H5Dwrite(sd_data_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dY->data);
  status &= H5Dclose(sd_data_id);
  status &= H5LTset_attribute_int(sd_group_id,H5_data_name,"index",&index, 1);
  status &= H5LTset_attribute_int(sd_group_id,H5_data_name,"major",&major, 1); // major experiment index, as presented in SBtab file
  status &= H5LTset_attribute_int(sd_group_id,H5_data_name,"minor",&minor, 1); // minor
  
  status &= H5Gclose(data_group_id);
  status &= H5Gclose(sd_group_id);
  assert(status==0);
  int current_index=index;
  if (status<0){
    fprintf(stderr,"[write_data_to_hdf5] something went wrong; overall H5 status=%i for simulation unit %i (Experiment %i.%i)\n",status,current_index,major,minor);
  }else{
    index++;
  }  
  return status;
}


gchar* get_reference(gchar *LinkKey, regmatch_t *match){
  int l=match->rm_eo - match->rm_so;
  int L=strlen(LinkKey);
  int s_offset=match->rm_so;
  gchar *s;
  //printf("[get_reference] linked key: «%s».\n",LinkKey);
  //fflush(stdout);
  if (l>0 && s_offset>0 && l <= L - s_offset){
    s=g_strndup(LinkKey+s_offset, l);
    //printf("[get_reference] «%s».\n",s);
    fflush(stdout);
  }else{
    //fprintf(stderr,"[get_reference] key: %s; match %i:%i (%i).\n",LinkKey,match->rm_so,match->rm_eo,L);
    //fprintf(stderr,"[get_reference] s_offset: %i, length: %i\n",s_offset,l);
    s=g_strdup("");
  }
  return s; 
}

gchar* dup_match(regmatch_t *match, char *source){
  regoff_t a,b;
  int length;
  gchar *s;
  a=match->rm_so;
  b=match->rm_eo;
  length=b-a;
  if (length>0)
    s=g_strndup(source+a,length);
  else
    s=NULL;
  return s;
}


int copy_match(char *sptr, regmatch_t *match, char *source, ssize_t N){
  regoff_t a,b;
  int length;
  a=match->rm_so;
  b=match->rm_eo;
  length=b-a;
  memset(sptr,'\0',length+1);
  if (N>length+1 && sptr!=NULL){
    strncpy(sptr,source+a,length);
    return EXIT_SUCCESS;
  } else {
    return EXIT_FAILURE;
  }  
}

gsl_matrix** get_data_matrix(sbtab_t *DataTable, sbtab_t *L1_OUT, int lflag){
  int i,i_r, i_c;
  GPtrArray *c,*dc;
  double val,dval=INFINITY;
  //double *DefaultData, *DefaultError;
  
  //printf("[get_data_matrix] checking whether DataTable exists.\n"); fflush(stdout);
  assert(DataTable!=NULL && L1_OUT!=NULL);
  guint F=L1_OUT->column[0]->len;
  guint T=DataTable->column[0]->len;
  gsl_matrix **Y;
  gchar *y,*dy,*s, *ds;
  Y=malloc(sizeof(gsl_matrix*)*2);
  
  for (i=0;i<2;i++) {
    printf("[get_data_matrix] allocating a %i×%i matrix.\n",T,F);
    Y[i]=gsl_matrix_alloc(T,F);
    assert(Y[i]!=NULL);
  }
  printf("[get_data_matrix] setting data matrix to default values (1 ± ∞).\n"); fflush(stdout);
  gsl_matrix_set_all(Y[DATA],1.0);
  if (lflag){
    gsl_matrix_set_all(Y[STDV],INFINITY);
  }else{ // so, lflag is FALSE, that means that it's a normalisation experiment; we avoid devision by infinity just to be sure:
    gsl_matrix_set_all(Y[STDV],1024);
  }


  // find the name of each outputs error quantifier:  
  GPtrArray *NoiseName;
  printf("[get_data_matrix] Found experiment %s with %i measurements of %i items.\n",DataTable->TableName,T,F);  fflush(stdout);
  for (i_c=0; i_c<F; i_c++){
    y=(gchar *) g_ptr_array_index(L1_OUT->column[0],i_c);
    NoiseName=sbtab_get_column(L1_OUT,"!ErrorName");
    if (NoiseName!=NULL){
      dy=g_ptr_array_index(NoiseName,i_c);
    }else{
      dy=g_strconcat("SD",y,NULL);
    }
    
    //printf("[get_data_matrix] reading the column %s ± %s\n",y,dy);  fflush(stdout);
    c=sbtab_get_column(DataTable,y);
    dc=sbtab_get_column(DataTable,dy);
    if (c!=NULL){
      for (i_r=0; i_r<T; i_r++){
	s = (gchar*) g_ptr_array_index(c,i_r);
	if (s!=NULL){
	  //printf("[get_data_matrix] row %i; got string «%s» ",i_r,s);
	  val=strtod(s,NULL);
	  gsl_matrix_set(Y[DATA],i_r,i_c,val);
	} 
	if (dc!=NULL) {
	  ds = (gchar*) g_ptr_array_index(dc,i_r);
	  assert(ds!=NULL);
	  dval=strtod(ds,NULL);
	  gsl_matrix_set(Y[STDV],i_r,i_c,dval);
	}
	//printf("..so %g ± %g\n",val,dval);
      }
    }
  }  
  return Y;  
}

gsl_vector* get_time_vector(sbtab_t *DataTable){
  int i;
  GPtrArray *t;
  gchar *s;
  double val;
  guint T;
  gsl_vector *time;
  //printf("[get_data_matrix] checking whether DataTable exists.\n"); fflush(stdout);
  assert(DataTable!=NULL);
  printf("[get_time_vector] getting measurement times.\n"); fflush(stdout);
  t=sbtab_get_column(DataTable,"!Time");
  if (t!=NULL){
    T=t->len;    
    time=gsl_vector_alloc(T);      
    printf("[get_time_vector] Found time column.\n");  fflush(stdout);
    for (i=0;i<T;i++){
      s = (gchar*) g_ptr_array_index(t,i);
      if (s!=NULL){
	val=strtod(s,NULL);
	gsl_vector_set(time,i,val);
      }
    }
  }else{
    time=NULL;
  }  
  return time;  
}

gsl_matrix* get_input_matrix(sbtab_t *DataTable, GPtrArray *input_ID, gsl_vector *default_input){
  int i,j;
  //printf("[get_data_matrix] checking whether DataTable exists.\n"); fflush(stdout);
  assert(DataTable!=NULL);
  guint N=DataTable->column[0]->len;
  guint nU=default_input->size;
  gchar *s=NULL, *ui=NULL;
  GPtrArray *u=NULL;
  gsl_matrix *U;
  gsl_vector_view u_row;
  double val;
  U=gsl_matrix_alloc(N,nU);
  for (i=0;i<N;i++){
    // copy all default values into matrix, then override later.
    u_row=gsl_matrix_row(U,i);
    gsl_vector_memcpy(&(u_row.vector),default_input);
  }
  assert(input_ID!=NULL);
  printf("[get_input_matrix] getting dose values in a Dose Response experiment %s.\n",DataTable->TableName); fflush(stdout);
  for (j=0;j<nU;j++){
    s=g_ptr_array_index(input_ID,j);
    if (s) u=sbtab_get_column(DataTable,s);
    if (u){
      printf("[get_input_matrix] Found input column %i (%s).\n",j,s);  fflush(stdout);
      for (i=0;i<N;i++){
	ui = (gchar*) g_ptr_array_index(u,i);
	if (ui){
	  // override defaults
	  val=strtod(ui,NULL);
	  gsl_matrix_set(U,i,j,val);
	}
      }
    }
  }
  return U;  
}


sbtab_t* parse_sb_tab(char *sb_tab_file){
  FILE* fid;
  int i,i_sep;
  size_t L=0;
  gchar *fs;
  //  char *key;
  char *s; // string to hold read in lines
  size_t n_s=DEFAULT_STR_LENGTH;
  ssize_t m_s;
  regex_t SBtab;
  regex_t RE_TableName, RE_TableTitle, RE_TableType, SBcomment, SBkeys, SBkey, SBlink;
  regex_t EmptyLine;
  gchar *TableName, *TableTitle, *TableType;
  regmatch_t match[4];
  regoff_t a;
  int r_status=0;
  gchar **keys;
  gchar *stem, *leaf;
  sbtab_t *sbtab;
  
  fs=g_strdup(";\t,");
  sbtab=NULL;
  s=malloc(sizeof(char)*n_s); // a buffer to read strings from file via getline
  r_status&=regcomp(&EmptyLine,"^[[:blank:]]*$",REG_EXTENDED);
  r_status&=regcomp(&SBcomment,"[%#]",REG_EXTENDED);
  r_status&=regcomp(&SBtab,"!!SBtab",REG_EXTENDED);  
  r_status&=regcomp(&RE_TableName,"TableName[[:blank:]]*=[[:blank:]]*'([^']+)'",REG_EXTENDED);
  r_status&=regcomp(&RE_TableTitle,"TableTitle[[:blank:]]*=[[:blank:]]*'([^']+)'",REG_EXTENDED);
  r_status&=regcomp(&RE_TableType,"TableType[[:blank:]]*=[[:blank:]]*'([^']+)'",REG_EXTENDED);
  r_status&=regcomp(&SBkeys,"![^!][[:alpha:]]",REG_EXTENDED);
  r_status&=regcomp(&SBkey,"(![[:alpha:]][[:alnum:]]*|>([[:alpha:]][[:word:]]*:)*([[:alpha:]][[:word:]])*)",REG_EXTENDED);
  r_status&=regcomp(&SBlink,">([[:alpha:]][[:alpha:][:digit:]]*:)*([[:alpha:]][[:alpha:][:digit:]]*)",REG_EXTENDED);
  
  if (r_status==0){
    printf("[parse_sb_tab] regular expressions created (EC=%i).\n",r_status);
  } else {
    fprintf(stderr,"[regcomp] returned %i\n",r_status);
    perror("[regcomp]");
  }
  
  fid=fopen(sb_tab_file,"r");
  if (fid==NULL){
    perror("[parse_sb_tab] ");
  }else{
    while (!feof(fid)){
      m_s = getline(&s,&n_s,fid);
      //printf("[parse_sb_tab] %i characters read.\n",m_s);
      if (m_s>0){
	if (regexec(&SBcomment,s,1,match,0)==0){
	  // remove comment from line
	  a=match[0].rm_so;
	  s[a]='\0';
	}
	if (regexec(&SBtab,s,0,NULL,0)==0){
	  i_sep=strcspn(s,fs);
	  L=strlen(s);
	  if (i_sep<L){
	    g_free(fs);
	    fs=g_strndup(&s[i_sep],1); //field separator
	  }
	  printf("[parse_sb_tab] separator: «%s».\n",fs);
	  if (regexec(&RE_TableName,s,2,match,0)==0){
	    TableName=g_strndup(s+match[1].rm_so,match[1].rm_eo-match[1].rm_so);
	    printf("TableName: «%s»\n",TableName); fflush(stdout);
	  } else {
	    fprintf(stderr,"error: TableName is missing.\n");
	    exit(-1);
	  }
	  if (regexec(&RE_TableType,s,2,match,0)==0){
	    TableType=g_strndup(s+match[1].rm_so,match[1].rm_eo-match[1].rm_so);
	    printf("TableType: «%s»\n",TableType); fflush(stdout);
	  }else {
	    fprintf(stderr,"warning: TableType is missing.\n");
	  }
	  if (regexec(&RE_TableTitle,s,2,match,0)==0){
	    TableTitle=g_strndup(s+match[1].rm_so,match[1].rm_eo-match[1].rm_so);
	    printf("TableTitle: «%s»\n",TableTitle); fflush(stdout);
	  }else {
	    fprintf(stderr,"warning: TableTitle is missing.\n");
	  }	
	} else if (regexec(&SBkeys,s,2,match,0)==0){
	  keys=g_strsplit_set(s,fs,-1);
	  int k=g_strv_length(keys);
	  //print all headers
	  printf("[parse_sb_tab] %i headers:",k); fflush(stdout);
	  for (i=0;i<k;i++) {
	    g_strstrip(keys[i]);
	    if (keys[i]!=NULL) printf("«%s» ",keys[i]);
	    else fprintf(stderr,"key[%i/%i] is NULL. ",i,k); 
	  }
	  printf("done.\n");
	  // check headers for links and save only the id in the link as hash.
	  for (i=0;i<k;i++) {
	    //printf("[parse_sb_tab] checking key %i of %i (%s); ",i,k,keys[i]);
	    r_status=regexec(&SBlink,keys[i],4,match,0);
	    fflush(stdout);
	    fflush(stderr);
	    if (r_status==0){
	      //printf("key[%i] is linked.\n",i);
	      fflush(stdout);
	      stem=get_reference(keys[i], &match[1]);
	      leaf=get_reference(keys[i], &match[2]);
	      if (leaf!=NULL){
		printf("\t[keys] link to table «%s», ID=«%s» found. I will use ID in hash table.\n",stem,leaf);
		g_free(keys[i]);
		keys[i]=leaf;
	      } else {
		fprintf(stderr,"[parse_sb_tab] could not dereference link. Will use the string as is.\n");
	      }
	    }else{
	      //printf("key[%i] is not linked.\n",i);
	    }	  
	  }
	  sbtab=sbtab_alloc(keys);
	} else if (regexec(&EmptyLine,s,0,NULL,0)==0){
	  printf("[parse_sb_tab] skipping empty line «%s».\n",s);
	} else {
	  assert(sbtab!=NULL);
	  //printf("[parse_sb_tab] append «%s» (m_s=%i)\n",s,m_s);
	  sbtab_append_row(sbtab,s,fs);
	}
      }
    }
    if (sbtab!=NULL){
      sbtab->TableTitle=TableTitle;
      sbtab->TableName=TableName;
      sbtab->TableType=TableType;
    }
  }  
  return sbtab;
}
