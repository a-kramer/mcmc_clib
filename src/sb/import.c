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
#include "data.h"
#include "sbtab.h"
#include "re.h"
#include "../app/event.h"
#include "../mcmc/model_parameters_smmala.h"
#include "../mcmc/ptype.h"
#include "event_glib.h"
/* very few strings have a fixed length: */

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


typedef struct {
  GPtrArray *id;
  union {
    gsl_vector *value;
    sbtab_t *table;
  };
} KeyValue_t;

gsl_matrix** get_data_matrix(sbtab_t *DataTable, sbtab_t *Output, int IsNormalisingExperiment);
gsl_matrix* get_input_matrix(sbtab_t *DataTable, GPtrArray *input_ID, gsl_vector *default_input);
herr_t write_prior_to_hdf5(hid_t file_id, gsl_vector *mu, void *S, int type);
int process_data_tables(hid_t file_id,  GPtrArray *sbtab,  GHashTable *sbtab_hash);
herr_t process_prior(hid_t file_id, GPtrArray *sbtab, GHashTable *sbtab_hash);

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
  regcomp(sb_tab_file,".*[.]tsv$",REG_EXTENDED);
  regcomp(h5_file,".*[.]h5$",REG_EXTENDED);
  sbtab=g_ptr_array_new_full(3,sbtab_free);
  sbtab_hash=g_hash_table_new(g_str_hash, g_str_equal);
  for (i=1;i<argc;i++){
    if (regexec(sb_tab_file,argv[i],0,NULL,0)==0) {
      printf("[%s] parsing sbtab file %s (arg %i of %i).\n",__func__,argv[i],i,argc);
      fflush(stdout);
      table=sbtab_from_tsv(argv[i]);
      if (table){
	g_ptr_array_add(sbtab,table);
	g_hash_table_insert(sbtab_hash,table->TableName,table);
	printf("[%s] added table «%s» with title: «%s».\n",__func__,table->TableName,table->TableTitle);
	fflush(stdout);
      }
    } else if (regexec(h5_file,argv[i],0,NULL,0)==0) {
      H5FileName=g_strdup(argv[i]);
    } else printf("unknown option «%s».\n",argv[i]);
  }
  if (H5FileName==NULL){
      H5FileName=g_strdup("ExperimentalData.h5");
  }
  guint n_tab=sbtab->len;
  printf("[%s] %i tables read: ",__func__,n_tab);
  fflush(stdout);
  for (i=0;i<n_tab;i++){
    table=g_ptr_array_index(sbtab,i);
    printf(" %s ",table->TableName);
    fflush(stdout);
  }
  printf("\n");
  fflush(stdout);
  hid_t file_id;
  file_id = H5Fcreate(H5FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  status |= process_data_tables(file_id, sbtab, sbtab_hash);
  status |= process_prior(file_id, sbtab, sbtab_hash);
  status |= H5Fclose(file_id);
  //process_input_table(H5FileName, sbtab, sbtab_hash);
  //process_prior(H5FileName, sbtab, sbtab_hash);
  return status;
}

gsl_vector* get_default_input(sbtab_t *Input){
  int j;
  double val;
  guint nU=table_length(Input);
  gsl_vector *default_input;
  default_input=gsl_vector_alloc(nU);
  gsl_vector_set_zero(default_input);
  GPtrArray *u_col;
  u_col=sbtab_find_column(Input,"!DefaultValue",NULL);
  printf("[%s] ",__func__);
  if (u_col){
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

gsl_vector *get_default_time(sbtab_t *E_table){
  gsl_vector *default_time;
  GPtrArray *default_time_column;
  default_time_column=sbtab_find_column(E_table,"!Time",NULL);
  int nE=table_length(E_table);
  int j;
  gchar *s;
  double t;
  default_time=gsl_vector_alloc(nE);
  printf("[%s] ",__func__);
  if (default_time_column){
    for (j=0;j<nE;j++) {
      s=g_ptr_array_index(default_time_column,j);
      t=strtod(s,NULL);
      gsl_vector_set(default_time,j,t);
      printf(" %g",t);
    }
    printf(";\n");
  } else {    
    gsl_vector_set_zero(default_time);
    printf("not defined, setting to zero.\n");
  }
  return default_time;
}

experiment_type* get_experiment_type(GPtrArray *E_type){
  int j;
  experiment_type *ExperimentType;
  gchar *Type;
  regex_t DoseResponse;
  regcomp(&DoseResponse,"Dose[[:blank:]]*Response", REG_EXTENDED | REG_ICASE);
  regex_t TimeSeries;
  regcomp(&TimeSeries,"Time[[:blank:]]*Series", REG_EXTENDED | REG_ICASE);
  int nE=E_type->len;
  printf("[%s] %i Experiments:",__func__,nE);
  ExperimentType=malloc(nE*sizeof(experiment_type));  
  for (j=0;j<nE;j++){
    Type=NULL;
    if (E_type){
      Type=g_ptr_array_index(E_type,j);
      if (Type){
	printf( "«%s»",Type);
	if (regexec(&DoseResponse, Type, 0, NULL, 0)==0){
	  ExperimentType[j]=dose_response;
	} else if (regexec(&TimeSeries, Type, 0, NULL, 0)==0){
	  ExperimentType[j]=time_series;  
	} else {
	  //fprintf(stderr,"[%s] This type is not known: «%s».\n",__func__,Type);
	  ExperimentType[j]=unknown_type;
	}
      } else {
	ExperimentType[j]=unknown_type;
	fprintf(stderr,"[%s] Experiment Type is set, but is not one of: «Time Series», «Dose Response»).\n",__func__);	
      }      
    }else{
      ExperimentType[j]=time_series;
      printf("[%s] Experiment Type is not defined, assuming time series data. Use !Type column in Table of Experiments (Possible Values: «Time Series», «Dose Response»).\n",__func__);
    }
    printf("(%i) ",ExperimentType[j]);
  }
  printf(";\n");
  regfree(&DoseResponse);
  regfree(&TimeSeries);
  fflush(stdout);
  return ExperimentType;
}


void override_default_input(gpointer data, gpointer user_data){
  KeyValue_t *input=user_data;
  sbtab_t *InputSet=data;
  GPtrArray *id=input->id;
  guint n=id->len;
  gint i,j;
  gint ri;
  gchar *ID;
  double value;
  GPtrArray *Value=sbtab_find_column(InputSet,"!DefaultValue !Value",NULL);
  assert(Value);
  assert(InputSet);
  for (i=0;i<n;i++){
    ID=g_ptr_array_index(id,i);
    ri=sbtab_find_row(InputSet,ID);
    if (ri>=0){
      value=strtod(g_ptr_array_index(Value,ri),NULL);
      gsl_vector_set(input->value,i,value);
    }
  }
}

GPtrArray* get_experiment_specific_inputs(sbtab_t *ExperimentTable, sbtab_t *Input, gsl_vector *default_input, GHashTable *sbtab_hash){
  assert(default_input);
  int j,k;
  gchar *s, *uid;
  GPtrArray *u_col;
  double val;
  GPtrArray *input_ID;
  GPtrArray *InputSetNames=sbtab_find_column(ExperimentTable,"!InputSet",NULL);
  
  int nU=default_input->size;
  int nE=table_length(ExperimentTable);
  GPtrArray *E_default_input=g_ptr_array_new_full(nE,gsl_vector_free);
  gsl_vector *u;
  
  gchar *InputSets;
  GPtrArray *InputTable;
  KeyValue_t input;
  printf("[%s] input has size %i.\n",__func__,nU);
  for (j=0;j<nE;j++) { // j is the major experiment index
    u=gsl_vector_alloc(nU);

    gsl_vector_memcpy(u,default_input);
    input_ID=sbtab_find_column(Input,"!ID",NULL);
    
    assert(input_ID);
    if (InputSetNames && j<InputSetNames->len){
      InputSets=(gchar*) g_ptr_array_index(InputSetNames,j);
      InputTable=InputSets?sbtab_get_tables(sbtab_hash,InputSets):NULL;
      if(InputTable){
	input.id=input_ID;
	input.value=u;
	g_ptr_array_foreach(InputTable,override_default_input,&input);
      }
    }
    for (k=0;k<nU;k++){
      uid=g_ptr_array_index(input_ID,k);
      u_col=sbtab_find_column(ExperimentTable,uid,">");
      if (u_col){
	s=g_ptr_array_index(u_col,j);
	val=strtod(s,NULL);
	gsl_vector_set(u,k,val);
      }
    }
    g_ptr_array_add(E_default_input,u);
  }
  return E_default_input;
}

sbtab_t* get_data_table(GHashTable *sbtab_hash, gchar *experiment_id, gchar *experiment_name){
  sbtab_t *DataTable;
  DataTable=sbtab_find(sbtab_hash,experiment_id);
  if (DataTable==NULL) {
    fprintf(stderr,"[%s] Could not find Measurement Table by ID:\t«TableName='%s'», \n\t\twill try searching by Name instead:\t«TableName='%s'».\n",__func__,experiment_id,experiment_name);
    if (experiment_name) DataTable=sbtab_find(sbtab_hash,experiment_name);
  }
  assert(DataTable);
  printf("[%s] found table: «%s» (TableName).\n",__func__,DataTable->TableName);
  return DataTable;
}

/* loads each data table in mentioned in the «Experiments» table */
GPtrArray* /* array of DataTables as referenced by TableName in the table of experiments */
get_data(GHashTable *sbtab_hash, /* hash table of sbtab tables */
	 sbtab_t *Experiments, /* the full experiment table */
	 int *lflag)/*an additional output: the likelihood contribution flag for each experiment*/{ 
  regex_t TrueRE, FalseRE;
  regcomp(&TrueRE,"true|1",REG_EXTENDED|REG_ICASE|REG_NOSUB);
  regcomp(&FalseRE,"false|0",REG_EXTENDED|REG_ICASE|REG_NOSUB);
  gchar *experiment_name, *experiment_id;
  gchar *flag;
  GPtrArray *ID=sbtab_find_column(Experiments,"!ID",NULL);
  GPtrArray *Name=sbtab_find_column(Experiments,"!Name",NULL);
  guint j,nE=ID->len;
  GPtrArray *LikelihoodFlag=sbtab_find_column(Experiments,"!Likelihood",NULL);
  GPtrArray *DataTable=g_ptr_array_new_full(nE,sbtab_free);
  assert(lflag);
  for (j=0;j<nE;j++) { // j is the major experiment index
    if (ID) experiment_id=(gchar*) g_ptr_array_index(ID,j);
    else experiment_id=(gchar*) g_ptr_array_index(Experiments->column[0],j);
    if (Name) {
      experiment_name=(gchar*) g_ptr_array_index(Name,j);
    } else {
      experiment_name=g_strdup(experiment_id);
    }
    printf("[%s] processing %s\n",__func__,experiment_id);
    g_ptr_array_add(DataTable,get_data_table(sbtab_hash, experiment_id, experiment_name));
    if (LikelihoodFlag){
      flag=g_ptr_array_index(LikelihoodFlag,j);
      lflag[j]=(regexec(&FalseRE,flag,0,NULL,0)!=0);
      if (!lflag[j]){
	printf("[%s] Experiment[%i] %s «%s» will be used for normalisation purposes only; it will not explicitly appear in the likelihood.\n",__func__,j,experiment_id,experiment_name);
      }
    } else {
      lflag[j]=1;
    }
  }  
  regfree(&FalseRE);
  regfree(&TrueRE);
  return DataTable;
}


void determine_target_and_operation(gpointer key, gpointer value, gpointer buffer){
  gchar *Key=key;
  GPtrArray *TargetName;
  gchar *Name;
  GPtrArray *C=value;
  gchar *c;
  guint N=C->len;
  int EC;
  EventEnvironment *SBEE=buffer;
  sbtab_t *Compound=sbtab_find(SBEE->sbtab_hash,"Compound Compounds");
  sbtab_t *Input=sbtab_find(SBEE->sbtab_hash,"Input Inputs");
  regex_t OP_Target_RE;
  regmatch_t m[3];
  gchar *OP,*Target;
  gint *op;
  effect_t effect;
  char errbuf[80];
  double val;
  int r=-1,i,j,k;
  GPtrArray *ID;
  gchar *id;
  op_t type=event_set; /* default action*/
  if (N>0){
    EC=regcomp(&OP_Target_RE,">([[:alpha:]]+):([[:alnum:]]+)",REG_EXTENDED);
    if (EC){
      regerror(EC, &OP_Target_RE, errbuf, 80); abort();
    }
    if (regexec(&OP_Target_RE, Key,3,m,0)==0){
      OP=dup_match(&m[1],Key);
      Target=dup_match(&m[2],Key);
      r=sbtab_get_row_index(Compound,Target);
      effect=event_affects_state;
      if (r<0){
	fprintf(stderr,"[%s] did not find «%s» in «%s»; all keys: \n",__func__,Target,Compound->TableName);
	ID=sbtab_find_column(Compound,"!ID",NULL);
	for (j=0;j<ID->len;j++){
	  id=g_ptr_array_index(ID,j);
	  fprintf(stderr,"\t%i «%s»\n",j,id);
	}
	r=sbtab_get_row_index(Input,Target);
	assert(r>=0);
	effect=event_affects_input;
      }
      g_array_append_val(SBEE->event->effect,effect);

      switch (effect){
      case event_affects_state:
	TargetName=sbtab_find_column(Compound,"!Name",NULL);
	break;
      case event_affects_input:
	TargetName=sbtab_find_column(Input,"!Name",NULL);
	break;
      default : fprintf(stderr,"[%s] effect unknown: %i.\n",__func__,effect); abort();	
      }
      Name=g_ptr_array_index(TargetName,r);
      if (Name){
	g_ptr_array_add(SBEE->event->gp_target_name,Name);
      }
      op=g_hash_table_lookup(SBEE->Operation,OP);
      type=op[0];
      g_array_append_val(SBEE->event->target,r);
      g_array_append_val(SBEE->event->Op,type);
      for (i=0;i<N;i++){
	c=g_ptr_array_index(C,i);
	if (c[0] == '>'){
	  /* a link */
	  k=sbtab_get_row_index(Input,&c[1]);
	  if (k<0){
	    fprintf(stderr,"[%s] event handling: could not find ID «%s» in Input table. A reference in that location «>ID» must refer to an input.\n",__func__,&c[1]);
	    abort();
	  } //else {
	    //fprintf(stdout,"[%s] event handling: found ID «%s» in Input table (row %i).\n",__func__,&c[1],k);
	  //}
	  assert(k>=0 && k < SBEE->Data->input->size);
	  val=gsl_vector_get(SBEE->Data->input,k);
	} else {
	  /* it's just a number then */
	  val=strtod(c,NULL);
	}
	//printf("[%s] value: %f.\n",__func__,val); fflush(stdout);
	g_array_append_val(SBEE->event->value,val);
      }
      g_free(OP);
      g_free(Target);
    } //else {
      //fprintf(stderr,"[%s] column title «%s» did not match the regular expression (>OPERATION:TARGETID).\n",__func__,Key); 
    //} 
    regfree(&OP_Target_RE);
  }
}
/* TODO rewrite the SBEE structure to contain substructures that can be allocated together consistently. a row can be represented well by glib structures, a numerical column also works well as a gsl object */
void event_interpret(gpointer event, gpointer buffer){
  EventEnvironment *SBEE=buffer;
  sbtab_t *sb_event=event;
  guint M=g_hash_table_size(sb_event->col);
  guint N=table_length(sb_event);
  guint L=M-1;
  guint i;
  guint k;
  GPtrArray *TimePoint=sbtab_find_column(sb_event,"!TimePoint",NULL);
  g_string_printf(SBEE->EventName,"%s_%s_%i",
		  sb_event->TableName,
		  SBEE->ExperimentName,
		  SBEE->ExperimentMinorIndex);
  if (TimePoint) L--;
  
  if (N>0 && L>0){
    SBEE->event=g_event_alloc(sbtab_column_to_gsl_vector(sb_event,"!Time"),M);
    g_hash_table_foreach(sb_event->col,determine_target_and_operation,SBEE);
    assert(SBEE->event->Op->len == SBEE->event->effect->len);
    assert(SBEE->event->gp_target_name->len == SBEE->event->Op->len);
    /* some output to check how the table was read:
     */
    size_t n=SBEE->event->target->len;
    size_t m=SBEE->event->time->size;
    
    //printf("[%s] in (%s) «%s» the (%li) targets are: ",__func__,SBEE->ExperimentName,SBEE->ExperimentID,n); fflush(stdout);
    for (i=0; i<n; i++){
      assert(SBEE->event->gp_target_name);
      assert(i<SBEE->event->gp_target_name->len);
      //printf("\t«%s» (row index: %i)\n",
      //	     (char*) g_ptr_array_index(SBEE->event->gp_target_name,i),
      //	     g_array_index(SBEE->event->target,int,i)); fflush(stdout);
    }
    // (1) make the values a gsl_matrix
    gsl_matrix *value_matrix=gsl_matrix_alloc(m,n);
    gsl_matrix_const_view g_mv_value=gsl_matrix_const_view_array((double*)(SBEE->event->value->data),n,m);
    // (2) transpose event matrix
    assert
      (gsl_matrix_transpose_memcpy
       (value_matrix,
	&(g_mv_value.matrix))==GSL_SUCCESS);
    // (3) write matrix to file
    //fprintf(stderr,"[%s] event «%s» value matrix: ",__func__,SBEE->EventName->str);
    //gsl_matrix_fprintf(stderr,value_matrix,"%g,");
    assert(SBEE->h5);
    SBEE->h5->size[0]=m;
    SBEE->h5->size[1]=n;
    printf("[%s] Experiment index: %i (%i.%i); Event Name: «%s»\n",__func__,SBEE->ExperimentIndex,SBEE->ExperimentMajorIndex,SBEE->ExperimentMinorIndex,SBEE->EventName->str); fflush(stdout);

    herr_t status;
    status=H5LTmake_dataset_double
      (SBEE->h5->group_id,
       SBEE->EventName->str,
       2, SBEE->h5->size,
       value_matrix->data);
    gsl_matrix_free(value_matrix);
    assert(status==0);
    status=H5LTset_attribute_int
      (SBEE->h5->group_id,
       SBEE->EventName->str,
       "AffectsMajorIndex",
       &(SBEE->ExperimentMajorIndex),1);
    assert(status==0);
    
    status=H5LTset_attribute_int
      (SBEE->h5->group_id,
       SBEE->EventName->str,
       "AffectsMinorIndex",
       &(SBEE->ExperimentMinorIndex),1);
    assert(status==0);
    
    status=H5LTset_attribute_int
      (SBEE->h5->group_id,
       SBEE->EventName->str,
       "index",
       &(SBEE->ExperimentIndex),1);
    assert(status==0);
    
    assert(SBEE->ExperimentName);
    status=H5LTset_attribute_string
      (SBEE->h5->group_id,
       SBEE->EventName->str,
       "ExperimentName",
       SBEE->ExperimentName);
    assert(status==0);
    
    status=H5LTset_attribute_string
      (SBEE->h5->group_id,
       SBEE->EventName->str,
       "ExperimentID",
       SBEE->ExperimentID);
    assert(status==0);

    assert(SBEE->event->Op->len>0);
    status=H5LTset_attribute_int
      (SBEE->h5->group_id,
       SBEE->EventName->str,
       "op",
       (int*) SBEE->event->Op->data,SBEE->event->Op->len);
    assert(status==0);

    assert(SBEE->event->effect->len);
    status=H5LTset_attribute_int
      (SBEE->h5->group_id,
       SBEE->EventName->str,
       "Effect",
       (int*) SBEE->event->effect->data, SBEE->event->effect->len);
    assert(status==0);

    status=H5LTset_attribute_double
      (SBEE->h5->group_id,
       SBEE->EventName->str,
       "Time",
       SBEE->event->time->data, SBEE->event->time->size);
    assert(status==0);
    
    /* Number of species can be different from SBtab, so matching by
       name is required */
    GString *TargetName = g_string_sized_new (20*n);
    for (i=0;i<n;i++) {
      g_string_append(TargetName,g_ptr_array_index(SBEE->event->gp_target_name,i));
      if (i<n-1) g_string_append_c(TargetName,' ');
    }
    //printf("[%s] TargetName: «%s»\n",__func__,TargetName->str); fflush(stdout);
    status=H5LTset_attribute_string
      (SBEE->h5->group_id,
       SBEE->EventName->str,
       "TargetName",
       TargetName->str);
    assert(status==0);
    g_string_free(TargetName,TRUE);

    GString *OP = g_string_sized_new(4*n);
    for (i=0;i<n;i++) {
      k=g_array_index(SBEE->event->Op,gint,i);
      g_string_append(OP,SBEE->OP_LABEL[k]);
      if (i<n-1) g_string_append_c(OP,' ');
    }
    // printf("[%s] Operation: «%s»\n",__func__,OP->str); fflush(stdout);
    status=H5LTset_attribute_string
      (SBEE->h5->group_id,
       SBEE->EventName->str,
       "Operation",
       OP->str);
    assert(status==0);
    g_string_free(OP,TRUE);
    GString *affected_table = g_string_sized_new(10*n);
    for (i=0;i<n;i++) {
      k=g_array_index(SBEE->event->effect,gint,i);
      switch(k){
      case event_affects_input:
	g_string_append(affected_table,"input");
	break;
      case event_affects_state:
	g_string_append(affected_table,"state");
	break;
      }
      if (i<n-1) g_string_append_c(affected_table,' ');
    }
    //printf("[%s] effect: «%s»\n",__func__,affected_table->str); fflush(stdout);
    status=H5LTset_attribute_string
      (SBEE->h5->group_id,
       SBEE->EventName->str,
       "affects",
       affected_table->str);
    assert(status==0);
    g_string_free(affected_table,TRUE);    
  } else {
    fprintf(stderr,"[%s] warning: empty event «%s»?\n",__func__,sb_event->TableTitle);
  }
  printf("[%s] done.\n",__func__); fflush(stdout);
}


/* this function gets each string valued cell in the !Event column and loads the events identified by that string.*/
void event_foreach_reference(gpointer PtrArrayElement, gpointer buffer){
  EventEnvironment *SBEE=buffer;
  GHashTable *sbtab_hash=SBEE->sbtab_hash;
  gchar *event_string=PtrArrayElement;
  GPtrArray *events=sbtab_get_tables(sbtab_hash,event_string);
  printf("[%s] processing experiment with %i events: «%s».\n",__func__,events->len,event_string);
  if (events->len>0){
    g_ptr_array_foreach(events,event_interpret,buffer);
  }
}

/**/
void event_foreach_experiment
(sbtab_t *Experiments,
 h5block_t *h5,
 GPtrArray *Data,
 sbtab_t *Input,
 map_t *IdxMap,
 GHashTable *sbtab_hash)
{
  GPtrArray *EventColumn=sbtab_find_column(Experiments,"!Event",NULL);
  gchar *e_str;
  sbtab_t *Parameter=sbtab_find(sbtab_hash,"Parameter Parameters");
  GPtrArray *ID,*Name;
  EventEnvironment SBEE;
  op_t OP[5]={event_set,
	      event_add,event_sub,
	      event_mul,event_div};
  gchar **OP_LABEL=g_strsplit("SET ADD SUB MUL DIV"," ",-1);
  int i,K=IdxMap->major->len;
  int major, minor;
  GString *EName=g_string_sized_new(40);
  if (EventColumn){
    assert(Parameter);
    SBEE.h5=h5;
    SBEE.sbtab_hash=sbtab_hash;
    ID=sbtab_find_column(Experiments,"!ID",NULL);
    Name=sbtab_find_column(Experiments,"!Name",NULL);
    SBEE.OP_LABEL=OP_LABEL;
    SBEE.Operation=g_hash_table_new(g_str_hash,g_str_equal);
    for (i=0;i<g_strv_length(OP_LABEL);i++){
      g_hash_table_insert(SBEE.Operation,OP_LABEL[i],&OP[i]);
    }
    
    for (i=0;i<K;i++){
      major=major_index(IdxMap,i);
      minor=minor_index(IdxMap,i);

      SBEE.ExperimentIndex=i;
      SBEE.ExperimentMajorIndex=major;
      SBEE.ExperimentMinorIndex=minor;
      SBEE.ExperimentName=g_ptr_array_index(Name,major);
      SBEE.ExperimentID=g_ptr_array_index(ID,major);
      SBEE.EventName=EName;
      SBEE.Data=g_ptr_array_index(Data,i);
      e_str=g_ptr_array_index(EventColumn,i);
      event_foreach_reference(e_str,&SBEE);
    }
  } else {
    printf("[%s] None of the experiments have events associated with them.\n",__func__);
  }
  g_string_free(EName,TRUE);
}

int process_data_tables(hid_t file_id,  GPtrArray *sbtab,  GHashTable *sbtab_hash){
  int status=EXIT_SUCCESS;
  sbtab_t *Input=sbtab_find(sbtab_hash,"Input input Inputs inputs");
  sbtab_t *Output=sbtab_find(sbtab_hash,"Output output L1OUT Outputs outputs");
  sbtab_t *Experiments=sbtab_find(sbtab_hash,"Experiments ExperimentTable TableOfExperiments");
  assert(Experiments && Input && Output);
  
  guint nE=table_length(Experiments);
  printf("[%s] %i Experiments.\n",__func__,nE);
  gsl_vector *default_time=get_default_time(Experiments);
  GPtrArray *EName=sbtab_find_column(Experiments,"!Name",NULL);
  GPtrArray *EType=sbtab_find_column(Experiments,"!Type",NULL);
  experiment_type *ExperimentType=get_experiment_type(EType);
  
  printf("[%s] found list of experiments: %s with %i entries\n",
	 __func__,Experiments->TableName,nE);


  gsl_vector *default_input=get_default_input(Input);
  GPtrArray *E_default_input=get_experiment_specific_inputs(Experiments, Input, default_input, sbtab_hash);

  int lflag[nE];
  GPtrArray *DataTable=get_data(sbtab_hash,Experiments,lflag);
  /* the above is as written inthe spreadsheet */
  printf("[%s] %i DataTables imported. Re-organizing into simulation packages.\n",__func__,DataTable->len); fflush(stdout);
  map_t *IdxMap=empty_map(nE);
  GPtrArray *D=unwrap_data
    (DataTable,
     lflag,
     ExperimentType,
     EName,
     E_default_input,
     default_time,
     Input,
     Output,
     IdxMap);
  printf("[%s] Data is re-packaged into %i simlation packages. Interpreting normalisation columns (!RelativeTo).\n",__func__,D->len); fflush(stdout);
  norm_t *N=normalisation(Experiments,ExperimentType,Output,DataTable,IdxMap);
  printf("[%s] normalisation read. Writing eerything to output file.\n",__func__); fflush(stdout);
  h5block_t *h5data=h5block_alloc(2);
  h5block_t *h5stdv=h5block_alloc(2);
  h5data->file_id=file_id;
  h5stdv->file_id=file_id;
  h5data->group_id=H5Gcreate2(file_id, "/data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5stdv->group_id=H5Gcreate2(file_id, "/stdv", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  data_attr_buffer_t buffer={h5data,h5stdv,IdxMap,N};

  g_ptr_array_foreach(D,write_data,&buffer);

  H5Gclose(h5data->group_id);
  H5Gclose(h5stdv->group_id);
  h5block_free(h5data);
  h5block_free(h5stdv);

  if (sbtab_find_column(Experiments,"!Event",NULL)){
    h5block_t *h5event=h5block_alloc(2);
    h5event->file_id=file_id;
    h5event->group_id=H5Gcreate2
      (file_id, "/event",
       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    event_foreach_experiment
      (Experiments, h5event,
       D, Input, IdxMap, sbtab_hash);
    H5Gclose(h5event->group_id);
    h5block_free(h5event);
  }
  map_free(IdxMap);
  return status;
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
  EC|=egrep_i(&LogType,"^(natural|base-e)?[[:blank:]]*(log(arithm)?)$|^ln$");
  EC|=egrep_i(&Log10Type,"^(decadic|base-10)[[:blank:]]*(log(arithm)?)$|^log10$");
  EC|=egrep_i(&LinType,"^lin(ear)?$");
  assert(EC==0);
 
  scale_type=gsl_vector_int_alloc((size_t) n);
  gsl_vector_int_set_all(scale_type,-1);
  for (i=0;i<n;i++){
    Type=g_ptr_array_index(Scale,i);
    if (regexec(&LogType, Type, 3, match, 0)==0){
      printf("[%s] (%i) «%s» matches RE for 'natural logarithm'\n",__func__,i,Type);
      gsl_vector_int_set(scale_type,i,LOG_SCALE);
    } else if (regexec(&Log10Type, Type, 3, match, 0)==0){
      printf("[%s] (%i) «%s» matches RE for 'base-10 logarithm'\n",__func__,i,Type);
      gsl_vector_int_set(scale_type,i,LOG10_SCALE);
    } else if (regexec(&LinType, Type, 3, match, 0)==0){
      printf("[%s] (%i) «%s» matches RE for 'linear scale'\n",__func__,i,Type);
      gsl_vector_int_set(scale_type,i,LIN_SCALE);
    } else {
      printf("[%s] This «!Scale[%i]» is unknown: «%s»\n",__func__,(int) i, Type);
      exit(-1);
    }
  }
  printf("[%s] ScaleType legend: natural logarithm = %i, base-10 logarithm = %i, linear = %i\n",__func__,LOG_SCALE,LOG10_SCALE,LIN_SCALE);
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
    /*   printf("[%s] (%i) prior value (parameter DefaultValues) seems to be in log-scale (natural logarithm); no scale conversion needed.\n",__func__,i); */
    }
  }
  return EXIT_SUCCESS;
}

herr_t process_prior(hid_t file_id, GPtrArray *sbtab, GHashTable *sbtab_hash){
  int i,n;
  printf("[%s] Looking up «Parameter» Table.\n",__func__);
  sbtab_t *Parameters;
  assert(sbtab_hash);
  Parameters=sbtab_find(sbtab_hash,"Parameter");
  assert(Parameters);
  GPtrArray *P_ID=sbtab_find_column(Parameters,"!ID",NULL);
  if (!P_ID) P_ID=sbtab_find_column(Parameters,"!Name",NULL);
  
  assert(P_ID);
  printf("[%s] Found an ID or Name column.\n",__func__);
  
  sbtab_t *Covariance, *Precision;
  Covariance=sbtab_find(sbtab_hash,"ParameterCovariance Covariance Sigma SIGMA cov");
  Precision=sbtab_find(sbtab_hash,"ParameterPrecision PriorPrecision Precision precision InvCovariance InverseCovariance inv(Sigma) inv(SIGMA)");
  printf("[%s] Checking Prior Type.\n",__func__);
  sbtab_t *cov[2]={Covariance, Precision};
  int any_cov_given=0;
  for (i=0;i<2;i++) any_cov_given+=(cov[i]!=NULL);
  printf("[%s] Any type of Covariance or Precision matrix given: %i.\n",__func__,any_cov_given); fflush(stdout);
  gchar *s;
  
  gsl_vector *mu=NULL;  
  gchar **ValueName;

  // get the scale (logarithmic or not)
  GPtrArray *Scale;
  Scale=sbtab_find_column(Parameters,"!Scale",NULL);
  gsl_vector_int *ScaleType; // 0 linear; 1 natural logarithm; 10 decadic logarithm;
  if (Scale){
    fprintf(stdout,"[%s] «!Scale» found in Parameter Table.\n",__func__);
    ScaleType=get_parameter_scale(Scale);
  }else{
    fprintf(stderr,"[%s] «!Scale» not found in Parameter Table, defaulting to LOG10 Scale.\n",__func__);
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
    printf("[%s] Using column «%s» for prior μ (mu).\n",__func__,ValueName[i]);
    g_strfreev(ValueName);
  }
  int D=mu->size;

  void *S; // S is either Precision, Covariance, or sigma;
  int type=PRIOR_IS_UNKNOWN;
  gsl_matrix *M=NULL;
  printf("[%s] prior has size %i.\n",__func__,D);

  if (any_cov_given){
    sbtab_t *C=NULL;
    type=PRIOR_IS_GAUSSIAN;
    ALSO(type,PRIOR_IS_MULTIVARIATE);

    M=gsl_matrix_alloc(D,D);
    gsl_matrix_set_identity(M);
    gsl_vector_view column;
    gsl_vector *c;
    if (Precision){
      ALSO(type,PRIOR_PRECISION_GIVEN);
      C=Precision;
    } else if (Covariance) {
      ALSO(type,PRIOR_SIGMA_GIVEN);
      C=Covariance;
    } else {
      printf("[%s] this prior specification {%i} is not handled yet.\n",__func__,any_cov_given);
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
      printf("[%s] using !Std column.\n",__func__);
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
	printf("[%s] Creating generalised Gaussian prior from !Min and !Max fields.\n",__func__);
	// apha is the width of the [min,max] interval 
	gsl_vector_add(alpha,par_max);
	gsl_vector_sub(alpha,par_min);
	gsl_vector_scale(alpha,0.5);
	S=alpha;
      }else{
	printf("[%s] neither !Std, nor (!Max,!Min) pairs were specified.\n",__func__);
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
  status |= H5LTmake_dataset_double(prior_group_id,"mu",1,&mu_size,mu->data);
  printf("[%s] prior type %i\n",__func__,type);
  if (PTYPE(type,PRIOR_IS_GAUSSIAN)){    
    if (PTYPE(type,PRIOR_IS_MULTIVARIATE)){
      printf("[%s] prior is multivariate Gaussian.\n",__func__); fflush(stdout);
      Sigma=(gsl_matrix*) S;      
      hsize_t Sigma_size[2];
      Sigma_size[0]=Sigma->size1;
      Sigma_size[1]=Sigma->size2;
      if (PTYPE(type,PRIOR_SIGMA_GIVEN))           Name=g_strdup("Sigma");
      else if (PTYPE(type,PRIOR_PRECISION_GIVEN))  Name=g_strdup("Precision");
      status |= H5LTmake_dataset_double(prior_group_id,Name,2,Sigma_size,Sigma->data);
    } else {
      printf("[write_prior_to_hdf5] prior is product of Gaussians.\n"); fflush(stdout);
      sigma=(gsl_vector*) S;
      if (PTYPE(type,PRIOR_IS_GENERALISED))        Name=g_strdup("alpha");
      else                                         Name=g_strdup("sigma");
      hsize_t sigma_size[1];
      sigma_size[0]=(hsize_t) sigma->size;
      status |= H5LTmake_dataset_double(prior_group_id,Name,1,sigma_size,sigma->data);
    }
    assert(prior_group_id>=0);
    printf("[%s] prior_group_id = %li.\n",__func__,prior_group_id);
    //status |= H5Gopen(file_id,"/prior",H5P_DEFAULT);
  } else {
    printf("[%s] Case not handled: %i\n",__func__,type);
    exit(-1);
  }
  if (Name) g_free(Name);
  status |= H5Gclose(prior_group_id);
  return status;
}

gsl_matrix** get_data_matrix(sbtab_t *DataTable, sbtab_t *Output, int lflag){
  int i,i_r, i_c;
  GPtrArray *c,*dc;
  double val,dval=INFINITY;
  assert(DataTable && Output);
  guint F=table_length(Output);
  guint T=table_length(DataTable);
  gsl_matrix **Y;
  gchar *y,*dy,*s, *ds;
  Y=malloc(sizeof(gsl_matrix*)*2);
  assert(F>0 && T>0);
  for (i=0;i<2;i++) {
    printf("[%s] allocating a %i×%i matrix.\n",__func__,T,F);
    Y[i]=gsl_matrix_alloc(T,F);
    assert(Y[i]);
  }
  printf("[%s] setting data matrix to default values (1 ± ∞).\n",__func__);
  fflush(stdout);
  gsl_matrix_set_all(Y[DATA],1.0);
  gsl_matrix_set_all(Y[STDV],INFINITY);
  GPtrArray *NoiseName;
  printf("[%s] Found experiment %s with %i measurements of %i items.\n",
	 __func__,DataTable->TableName,T,F);
  fflush(stdout);
  for (i_c=0; i_c<F; i_c++){
    y=(gchar *) g_ptr_array_index(Output->column[0],i_c);
    NoiseName=sbtab_find_column(Output,"!ErrorName",NULL);
    if (NoiseName){
      dy=g_ptr_array_index(NoiseName,i_c);
    }else{
      dy=g_strconcat("SD",y,NULL);
    }
    //printf("[%s] reading the column %s ± %s\n",__func__,y,dy);  fflush(stdout);
    c=sbtab_find_column(DataTable,y,">");
    dc=sbtab_find_column(DataTable,dy,NULL);
    if (c){
      for (i_r=0; i_r<T; i_r++){
	s = (gchar*) g_ptr_array_index(c,i_r);
	if (s){
	  printf("[%s] row %i; got string «%s» ",__func__,i_r,s);
	  val=strtod(s,NULL);
	  gsl_matrix_set(Y[DATA],i_r,i_c,val);
	} 
	if (dc) {
	  ds = (gchar*) g_ptr_array_index(dc,i_r);
	  assert(ds);
	  dval=strtod(ds,NULL);
	  gsl_matrix_set(Y[STDV],i_r,i_c,dval);
	}
	printf("..so %g ± %g\n",val,dval);
      }
    }
  }  
  return Y;  
}

gsl_matrix* get_input_matrix(sbtab_t *DataTable, GPtrArray *input_ID, gsl_vector *default_input){
  int i,j;
  //printf("[%s] checking whether DataTable exists.\n",__func__); fflush(stdout);
  assert(DataTable);
  guint N=DataTable->column[0]->len;
  guint nU=default_input->size;
  gchar *s=NULL, *ui=NULL, *uid_j;
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
  assert(input_ID);
  printf("[%s] getting dose values in a Dose Response experiment %s.\n",__func__,DataTable->TableName); fflush(stdout);
  for (j=0;j<nU;j++){
    uid_j = g_ptr_array_index(input_ID,j);
    u=sbtab_find_column(DataTable,uid_j,">");
    if (u){
      printf("[%s] Found input column %i (%s).\n",__func__,j,s);  fflush(stdout);
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

