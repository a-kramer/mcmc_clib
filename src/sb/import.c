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
#include "sbtab.h"
#include "../mcmc/model_parameters_smmala.h"

#define DEFAULT_STR_LENGTH 128

sbtab_t* parse_sb_tab(char *);
gsl_matrix** get_data_matrix(sbtab_t *DataTable,sbtab_t *L1_OUT);
int write_data_to_hdf5(char *file_name, gsl_matrix **Y, gsl_matrix **dY, int nE);
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

write_input_to_hdf5(gchar *H5FileName, gsl_matrix *u){
  int i;
  hid_t       file_id;   /* file identifier */
  herr_t      status;
  hid_t data_group_id, sd_group_id, dataspace_id, dataset_id, sd_data_id;  
  hsize_t size[2];
  char H5_data_name[24];
  status=EXIT_SUCCESS;
  file_id = H5Fopen(H5FileName, H5F_ACC_RDWR, H5P_DEFAULT);
  group_id = H5G(file_id, "/input", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /* do stuff */
  H5Gclose(group_id);
  H5Fclose(file_id);
}

int process_input_table(gchar *H5FileName,  GPtrArray *sbtab,  GHashTable *sbtab_hash){
  sbtab_t *E_table;
  guint nE;
  gchar *experiment_id,*s;
  GPtrArray *ID, *InputID, *ExperimentID, *E_U; // E_U is the experiment's input column
  sbtab_t *L1_OUT, *Input;
  gsl_matrix *u;
  int i,j,nE,nU; // number of Experiments and number of Iputs
  gchar *s,*;
  int status=EXIT_SUCCESS;
  double val;
  Input=g_hash_table_lookup(sbtab_hash,"Input");
  E_table=g_hash_table_lookup(sbtab_hash,"Experiments");
  if (E_table==NULL){
    E_table=g_hash_table_lookup(sbtab_hash,"TableOfExperiments");
  }
  if(E_table!=NULL && Input!=NULL){
    ID=E_table->column[0];
    nE=ID->len;
    ExperimentID=ID;

    ID=Input->column[0];
    nU=ID->len;
    InputID=ID;
    u=gsl_matrix_alloc(nE,nU);
    gsl_matrix_set_zero(u);
    for (i=0;i<nE;i++){      
      for (j=0;j<nU;j++){
	U=g_ptr_array_index(InputID,j);
	assert(U!=NULL);
	E_U=sbtab_get_column(E_table,U);
	s=g_ptr_array_index(E_U,i);
	val=strtod(s,NULL);
	gsl_matrix_set(u,i,j,val);
      }
    }
  }else{
    fprintf(stderr,"[process_input_table] Table missing (Experiments or Inputs)\n");
    exit(-1);
  }
  if (u!=NULL){
    write_input_to_hdf5(H5FileName,u);
  }
  return status;
}

int process_data_tables(gchar *H5FileName,  GPtrArray *sbtab,  GHashTable *sbtab_hash){
  // convert all data matrices of numbers into actual double/int numbers
  // 1. find all experiment quantity matrices
  sbtab_t *E_table;
  guint nE;
  gchar *experiment_id,*s;
  gsl_matrix **Y_dY;
  gsl_matrix **Y;
  gsl_matrix **dY;
  double val;
  sbtab_t *DataTable;
  sbtab_t *L1_OUT;
  GPtrArray *ID; 
  int status=EXIT_SUCCESS;
  
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
    nE=E_table->column[0]->len; // list of experiments -> number of experiments
    ID=g_hash_table_lookup(E_table->col,"!ID");
    if (ID!=E_table->column[0]) printf("!ID is not first column!\n");
    for (j=0;j<nE;j++){    
      experiment_id=(gchar*) g_ptr_array_index(E_table->column[0],j);
      printf("GPtrArray[%i]: %s\n",j,experiment_id);
      experiment_id=(gchar*) g_ptr_array_index(ID,j);
      printf("ID[%i]: %s\n",j,experiment_id);      
    }
    printf("[main] found list of experiments: %s with %i entries\n",E_table->TableName,nE);
    Y=malloc(sizeof(gsl_matrix*)*nE);
    dY=malloc(sizeof(gsl_matrix*)*nE);
    for (j=0;j<nE;j++) {
      experiment_id=(gchar*) g_ptr_array_index(E_table->column[0],j);
      printf("[main] processing %s\n",experiment_id);
      DataTable=g_hash_table_lookup(sbtab_hash,experiment_id);
      if (DataTable!=NULL && L1_OUT!=NULL){
	Y_dY=get_data_matrix(DataTable,L1_OUT);
	Y[j]=Y_dY[0];
	dY[j]=Y_dY[1];
      } else{
	fprintf(stderr,"Either L1 Output or DataTable %s missing (NULL).\n",experiment_id);
	status&=EXIT_FAILURE;
      }
    }
    //printf("result of reading matrices (Y,dY):\n");
    //for (j=0;j<nE;j++) gsl_matrix_fprintf(stdout,Y[j],"%g,");
  }
  if (Y!=NULL && dY!=NULL){
    write_data_to_hdf5(H5FileName,Y,dY,nE);
  }
  return status;
}


int write_data_to_hdf5(char *file_name, gsl_matrix **Y, gsl_matrix **dY, int nE){
  int i;
  hid_t       file_id;   /* file identifier */
  herr_t      status;
  hid_t data_group_id, sd_group_id, dataspace_id, dataset_id, sd_data_id;  
  hsize_t size[2];
  char H5_data_name[24];
  /* Create a new file using default properties. */
  status=EXIT_SUCCESS;
  file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  data_group_id = H5Gcreate(file_id, "/data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  sd_group_id = H5Gcreate(file_id, "/sd_data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for (i=0;i<nE;i++){
    size[0]=Y[i]->size1;
    size[1]=Y[i]->size2;
    dataspace_id = H5Screate_simple(2, size, NULL);
    // Y
    sprintf(H5_data_name,"data_block_%i",i);
    dataset_id = H5Dcreate2(data_group_id, H5_data_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status &= H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Y[i]->data);    
    status &= H5Dclose(dataset_id);
    // dY
    sprintf(H5_data_name,"sd_data_block_%i",i);
    sd_data_id = H5Dcreate2(sd_group_id, H5_data_name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status &= H5Dwrite(sd_data_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dY[i]->data);    
    status &= H5Dclose(sd_data_id);    
  }
  /* Terminate access to the file. */
  status &= H5Gclose (data_group_id);
  status &= H5Gclose (sd_group_id);
  status &= H5Fclose(file_id);
  if (status!=EXIT_SUCCESS){
    fprintf(stderr,"[write HDF5] something went wrong; overall status=%i\n",status);
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
  gsl_matrix_set_all(Y[0],1.0);
  gsl_matrix_set_all(Y[1],INFINITY);


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
	  gsl_matrix_set(Y[0],i_r,i_c,val);
	} 
	if (dc!=NULL) {
	  ds = (gchar*) g_ptr_array_index(dc,i_r);
	  assert(ds!=NULL);
	  dval=strtod(ds,NULL);
	  gsl_matrix_set(Y[1],i_r,i_c,dval);
	}
	//printf("..so %g ± %g\n",val,dval);
      }
    }
  }  
  return Y;  
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
