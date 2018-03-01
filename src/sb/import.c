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

#define DEFAULT_STR_LENGTH 128
#define UNKNOWN_TYPE 0
#define DOSE_RESPONSE 1
#define TIME_SERIES 2

#define DATA 0
#define STDV 1

sbtab_t* parse_sb_tab(char *);
gsl_matrix** get_data_matrix(sbtab_t *DataTable,sbtab_t *L1_OUT);
gsl_vector* get_time_vector(sbtab_t *DataTable);
gsl_matrix* get_input_matrix(sbtab_t *DataTable, GPtrArray *input_ID, gsl_vector *default_input);
int write_data_to_hdf5(hid_t file_id, gsl_matrix *Y, gsl_matrix *dY, gsl_vector *time, gsl_vector *input, int major, int minor);
int process_data_tables(gchar *H5FileName,  GPtrArray *sbtab,  GHashTable *sbtab_hash);

int main(int argc, char*argv[]){
  int i,j;
  int status=EXIT_SUCCESS;
  regex_t *sb_tab_file, *h5_file;
  sbtab_t *table;
  GPtrArray *sbtab;
  GHashTable *sbtab_hash;
  char *sptr;
  gchar *H5FileName=NULL;
  
  sptr=malloc(sizeof(char)*64);
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
  status&=process_data_tables(H5FileName, sbtab, sbtab_hash);
  //process_input_table(H5FileName, sbtab, sbtab_hash);
  //process_prior(H5FileName, sbtab, sbtab_hash);
  return status;
}

int process_data_tables(gchar *H5FileName,  GPtrArray *sbtab,  GHashTable *sbtab_hash){
  // convert all data matrices of numbers into actual double/int numbers
  // 1. find all experiment quantity matrices
  sbtab_t *E_table, *Input;  
  gchar *experiment_id, *s, *input_id, *type;
  double default_time=1.0; // arbitrary  
  gsl_vector *default_input;
  gsl_vector *E_default_input;
  gsl_vector *time;
  gsl_vector *input;
  gsl_vector **u; //inputs
  gsl_matrix **Y_dY; // to make one return pointer possible: Y_dY[0] is data, Y_dY[1] is standard deviation
  gsl_matrix **Y;    // vector of matrices, one per input line
  gsl_matrix **dY;
  double val;
  sbtab_t *DataTable;
  sbtab_t *L1_OUT;
  GPtrArray *ID, *E_type, *input_ID, *default_time_column; 
  int status=EXIT_SUCCESS;
  int i,j,k;
  hid_t file_id;
  // find Input Table
  Input=g_hash_table_lookup(sbtab_hash,"Input");
  assert(Input!=NULL);
  // find all output names
  L1_OUT=g_hash_table_lookup(sbtab_hash,"L1_OUT");
  if (L1_OUT==NULL){
    L1_OUT=g_hash_table_lookup(sbtab_hash,"Output");
  }
  assert(L1_OUT!=NULL);

  E_table=g_hash_table_lookup(sbtab_hash,"Experiments");
  if (E_table==NULL){
    E_table=g_hash_table_lookup(sbtab_hash,"TableOfExperiments");
  }
  
  if (E_table!=NULL){
    guint nE=E_table->column[0]->len; // list of experiments -> number of experiments
    ID=g_hash_table_lookup(E_table->col,"!ID");
    assert(ID!=NULL);
    if (ID!=E_table->column[0]) printf("!ID is not first column!\n");
    for (j=0;j<nE;j++){    
      experiment_id=(gchar*) g_ptr_array_index(E_table->column[0],j);
      printf("GPtrArray[%i]: %s\n",j,experiment_id);
      experiment_id=(gchar*) g_ptr_array_index(ID,j);
      printf("ID[%i]: %s\n",j,experiment_id);
      // default measurement time?
      default_time_column=sbtab_get_column(E_table,"!Time");
    }
    printf("[main] found list of experiments: %s with %i entries\n",E_table->TableName,nE);
    u=malloc(sizeof(gsl_vector*)*nE);
    guint nU=Input->column[0]->len;
    default_input=gsl_vector_alloc(nU);
    gsl_vector_set_zero(default_input);
    GPtrArray *u_col;
    u_col=sbtab_get_column(Input,"!DefaultValue");
    if (u_col!=NULL){
      for (j=0;j<nU;j++){
	val=strtod(g_ptr_array_index(u_col,j),NULL);
	gsl_vector_set(default_input,j,val);
      }
    }
    Y=malloc(sizeof(gsl_matrix*)*nE);
    dY=malloc(sizeof(gsl_matrix*)*nE);
    E_default_input=gsl_vector_alloc(nU);
    input=gsl_vector_alloc(nU);
    for (j=0;j<nE;j++) { // j is the major experiment index      
      experiment_id=(gchar*) g_ptr_array_index(E_table->column[0],j);
      if (default_time_column!=NULL) {
	default_time=strtod(g_ptr_array_index(default_time_column,j),NULL);
      } else {
	default_time=1;
      }
      printf("[main] processing %s\n",experiment_id);
      input_ID=sbtab_get_column(Input,"!ID");
      assert(input_ID!=NULL);
      gsl_vector_memcpy(E_default_input,default_input);
      for (k=0;k<nU;k++){
	input_id=g_ptr_array_index(input_ID,k);	
	u_col=sbtab_get_column(E_table,input_id);
	if (u_col!=NULL){
	  s=g_ptr_array_index(u_col,j);	
	  val=strtod(s,NULL);
	  gsl_vector_set(E_default_input,k,val);
	}
      }
      DataTable=g_hash_table_lookup(sbtab_hash,experiment_id);
      file_id = H5Fcreate(H5FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

      if (DataTable!=NULL && L1_OUT!=NULL){
	E_type=sbtab_get_column(E_table,"!Type");
	int ExperimentType=UNKNOWN_TYPE;
	if (E_type!=NULL){
	  regex_t DoseResponse;
	  regcomp(&DoseResponse,"[[:blank:]]*[Dd]ose[[:blank:]]*[Rr]esponse", REG_EXTENDED);
	  type=g_ptr_array_index(E_type,j);
	  if (type!=NULL && regexec(&DoseResponse, type, 0, NULL, 0)==0){
	    ExperimentType=DOSE_RESPONSE;
	  } else {
	    ExperimentType=TIME_SERIES;  
	  }
	}else{
	  // default case
	  ExperimentType=TIME_SERIES;
	}
	switch(ExperimentType){
	case UNKNOWN_TYPE:
	  fprintf(stderr,"[process_data_table] error: unknown experiment type\n");
	  exit(-1);
	  break;
	case TIME_SERIES:
	  Y_dY=get_data_matrix(DataTable,L1_OUT);
	  time=get_time_vector(DataTable);
	  write_data_to_hdf5(file_id,Y_dY[DATA],Y_dY[STDV],time,E_default_input,j,0);
	  break;
	case DOSE_RESPONSE:
	  Y_dY=get_data_matrix(DataTable,L1_OUT);
	  gsl_matrix *input_block;
	  input_block=get_input_matrix(DataTable,input_ID,E_default_input);
	  gsl_vector_view time_view;
	  time_view=gsl_vector_view_array(&default_time,1);
	  int N=Y_dY[DATA]->size1;
	  int M=Y_dY[DATA]->size2;
	  gsl_matrix_view data_sub, sd_data_sub;
	  gsl_vector_view input_row;
	  for (i=0;i<N;i++){ // i is the minor "simulation unit" index within a dose response experiment
	    data_sub=gsl_matrix_submatrix(Y_dY[DATA],i,0,1,M);
	    sd_data_sub=gsl_matrix_submatrix(Y_dY[STDV],i,0,1,M);
	    input_row=gsl_matrix_row(input_block,i);
	    write_data_to_hdf5(file_id,&(data_sub.matrix),&(sd_data_sub.matrix),&(time_view.vector),&(input_row.vector),j,i);
	  }
	  break;
	}
      } else{
	fprintf(stderr,"Either L1 Output or DataTable %s missing (NULL).\n",experiment_id);
	status&=EXIT_FAILURE;
      }
    }
    //printf("result of reading matrices (Y,dY):\n");
    //for (j=0;j<nE;j++) gsl_matrix_fprintf(stdout,Y[j],"%g,");
  }
  status &= H5Fclose(file_id);
  return status;
}

int write_data_to_hdf5(hid_t file_id, gsl_matrix *Y, gsl_matrix *dY, gsl_vector *time, gsl_vector *input, int major, int minor){
  int i;
  herr_t      status;
  hid_t data_group_id, sd_group_id, dataspace_id, dataset_id, sd_data_id;  
  hsize_t size[2];
  char H5_data_name[24];
  static int index=0;
  /* Create a new file using default properties. */
  status=EXIT_SUCCESS;
  data_group_id = H5Gcreate(file_id, "/data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  sd_group_id = H5Gcreate(file_id, "/sd_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  // write data and standard deviation to file
    size[0]=Y->size1;
    size[1]=Y->size2;
    dataspace_id = H5Screate_simple(2, size, NULL);
    // Y
    sprintf(H5_data_name,"data_block_%i",index);
    dataset_id = H5Dcreate2(data_group_id, H5_data_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status &= H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Y->data);
    status &= H5LTset_attribute_int(dataset_id,H5_data_name,"index", &i, 1);
    status &= H5Dclose(dataset_id);
    status &= H5LTset_attribute_int(dataset_id, H5_data_name, "major", &major, 1); // major experiment index, as presented in SBtab file
    status &= H5LTset_attribute_int(dataset_id, H5_data_name, "minor", &minor, 1); // minor experiment index, within a dose response experiment
    // dY
    sprintf(H5_data_name,"sd_data_block_%i",index);
    sd_data_id = H5Dcreate2(sd_group_id, H5_data_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status &= H5Dwrite(sd_data_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dY->data);    
    status &= H5Dclose(sd_data_id);    
    // time
    status &= H5LTset_attribute_double(dataset_id, H5_data_name, "time", time->data, time->size);
    // input
    status &= H5LTset_attribute_double(dataset_id, H5_data_name, "input", input->data, input->size);
  /* Terminate access to the file. */
  status &= H5Gclose (data_group_id);
  status &= H5Gclose (sd_group_id);

  if (status!=EXIT_SUCCESS){
    fprintf(stderr,"[write HDF5] something went wrong; overall H5 status=%i\n",status);
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

int copy_match(char *sptr, regmatch_t *match, char *source, ssize_t N){
  regoff_t a,b;
  int i;
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

gsl_matrix** get_data_matrix(sbtab_t *DataTable, sbtab_t *L1_OUT){
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
  gsl_matrix_set_all(Y[STDV],INFINITY);


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
  double default_time=1.0;
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
    printf("[get_time_vector] did not find a time column, using internal default measurement time: %f.\n",default_time);
    T=1;
    time=gsl_vector_alloc(T);
    gsl_vector_set(time,0,default_time);
  }  
  return time;  
}

gsl_matrix* get_input_matrix(sbtab_t *DataTable, GPtrArray *input_ID, gsl_vector *default_input){
  int i,j;
  //printf("[get_data_matrix] checking whether DataTable exists.\n"); fflush(stdout);
  assert(DataTable!=NULL);
  guint N=DataTable->column[0]->len;
  guint nU=default_input->size;
  gchar *s, *ui;
  GPtrArray *u;
  gsl_matrix *U;
  gsl_vector_view u_row;
  double val;
  U=gsl_matrix_alloc(N,nU);
  for (i=0;i<N;i++){
    // copy all default values into matrix, then override later.
    u_row=gsl_matrix_row(U,i);
    gsl_vector_memcpy(&(u_row.vector),default_input);
  }
  printf("[get_input_matrix] getting dose values in a Dose Response experiment %s.\n",DataTable->TableName); fflush(stdout);
  for (j=0;j<nU;j++){
    s=g_ptr_index(input_ID,j);
    u=sbtab_get_column(DataTable,s);
    if (u!=NULL){
      printf("[get_input_matrix] Found input column %i (%s).\n",j,s);  fflush(stdout);
      for (i=0;i<N;i++){
	ui = (gchar*) g_ptr_array_index(u,i);
	if (ui!=NULL){
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
  size_t n_s=DEFAULT_STR_LENGTH, m_s;
  regex_t SBtab;
  regex_t RE_TableName, RE_TableTitle, RE_TableType, SBcomment, SBkeys, SBkey, SBlink;
  regex_t EmptyLine;
  gchar *TableName, *TableTitle, *TableType;
  regmatch_t match[4];
  regoff_t a,b;
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
  r_status&=regcomp(&SBlink,">([[:alpha:]][[:alpha:][:digit:]]*:)*([[:alpha:]][[[:alpha:][:digit:]]]*)",REG_EXTENDED);
  
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
	  fprintf(stderr,"TableName is missing.\n");
	  exit(-1);
	}
	if (regexec(&RE_TableType,s,2,match,0)==0){
	  TableType=g_strndup(s+match[1].rm_so,match[1].rm_eo-match[1].rm_so);
	  printf("TableType: «%s»\n",TableType); fflush(stdout);
	}else {
	  fprintf(stderr,"TableType is missing.\n");
	  exit(-1);
	}
	if (regexec(&RE_TableTitle,s,2,match,0)==0){
	  TableTitle=g_strndup(s+match[1].rm_so,match[1].rm_eo-match[1].rm_so);
	  printf("TableTitle: «%s»\n",TableTitle); fflush(stdout);
	}else {
	  fprintf(stderr,"TableTitle is missing.\n");
	  exit(-1);
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
	sbtab_append_row(sbtab,s);
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
