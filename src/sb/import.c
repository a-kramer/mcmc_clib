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
#include "sbtab.h"
#include "../mcmc/model_parameters_smmala.h"

#define DEFAULT_STR_LENGTH 128

sbtab_t* parse_sb_tab(char *);
 
int main(int argc, char*argv[]){
  int i;
  regex_t *sb_tab_file;
  sbtab_t *table;
  GPtrArray *sbtab;
  GHashTable *sbtab_hash;
  char *sptr;
  
  sptr=malloc(sizeof(char)*64);
  sb_tab_file=malloc(sizeof(regex_t));
  regcomp(sb_tab_file,".*[.][tc]sv",REG_EXTENDED);
  sbtab=g_ptr_array_new_full(3,sbtab_free);
  sbtab_hash=g_hash_table_new(g_str_hash, g_str_equal);
  for (i=1;i<argc;i++){
    if (regexec(sb_tab_file,argv[i],0,NULL,0)==0) {
      table=parse_sb_tab(argv[i]);
      g_ptr_array_add(sbtab,table);
      g_hash_table_insert(sbtab_hash,table->TableName,table);
    }
    else printf("unknown option «%s».\n",argv[i]);
  }
  guint n=sbtab->len;
  printf("[main] %i tables read.\n",n);
  // convert all data matrices of numbers into actual double/int numbers
  // 1. find all experiment quantity matrices
  regex_t sb_ListOfExperiments;
  regmatch_t match[5];
  
  guint nE;
  guint T,N;
  gchar *id,*s;

  gsl_matrix **Y;
  gsl_matrix **dY;
  double val;
  sbtab_t *DataTable;
  sbtab_t *L1_OUT;
  
  
  // find all output names
  L1_OUT=g_hash_table_lookup(sbtab_hash,"L1_OUT");
  N=L1_OUT->column[0]->len;
  regcomp(&sb_ListOfExperiment,"(([Ll]ist|[Tt]able)Of)?[Ee]xp(eriments)?",REG_EXTENDED);
  
  for (i=0;i<n;i++){
    table=(sbtab_t*) g_ptr_array_index(sbtab,i);
    if (regexec(sb_ListOfExperiments, table->TableName, 5, match , 0)==0){
      printf("found list of experiments: %s\n",table->TableName);
      nE=table->column[0]->len; // list of experiments -> number of experiments
      Y=malloc(sizeof(gsl_matrix)*nE);
      dY=malloc(sizeof(gsl_matrix)*nE);
      for (j=0;j<nrows;j++) {
	id=(gchar*) g_ptr_array_index(table->column[0],j);
	DataTable=g_hash_table_lookup(sbtab_hash,id);
	get_data_matrix(DataTable,L1_OUT,Y[j],dY[j]);
    }
  }
return EXIT_SUCCESS;
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

int get_data_matrix(sbtab_t *DataTable,sbtab_t *L1_OUT, gsl_matrix *Y, gsl_matrix *dY){
  int status=EXIT_SUCCESS;
  int i_r, i_c;
  guint T,N;
  gchar *id;
  regex_t LinkedHeader;
  
  regcomp(&LinkedHeader,"^>(([[:alpha:]][[:word:]]*):)*([[:alpha:]][[:word:]]+)$",REG_EXTENDED);
  N=L1_OUT->column[0]->len;
  if (DataTable!=NULL){	  
    T=DataTable->column[0]->len;
    Y=gsl_matrix_alloc(T,N);	  
    for (i_r=0;i_r<T;i_r++){
      regexec();
      s=sbtab_get(DataTable,id,);
      val=strtod();
    }
    printf("Found experiment %s with %i measurements.\n",id,T);	  
  }else{
    fprintf(stderr,"Experiment %s could not be found.\n",id);
    status&=EXIT_FAILURE;
  }
  return status;  
}


sbtab_t* parse_sb_tab(char *sb_tab_file){
  FILE* fid;
  int i,i_sep;
  char fs=';';
  char *key;
  char *s; // string to hold read in lines
  size_t n_s=DEFAULT_STR_LENGTH, m_s, length=0;
  regex_t SBtab;
  regex_t RE_TableName, RE_TableTitle, RE_TableType, SBcomment, SBkeys, SBkey, SBlink;
  regex_t EmptyLine;
  gchar *TableName, *TableTitle, *TableType;
  regmatch_t match[2];
  regoff_t a,b;
  int r_status=0;
  gchar **keys;
  gchar *stem, *leaf;
  sbtab_t *sbtab;
  sbtab=NULL;
  s=malloc(sizeof(char)*n_s);
  r_status&=regcomp(&EmptyLine,"^[[:blank:]]*$",REG_EXTENDED);
  r_status&=regcomp(&SBcomment,"[%#]",REG_EXTENDED);
  r_status&=regcomp(&SBtab,"!!SBtab",REG_EXTENDED);  
  r_status&=regcomp(&RE_TableName,"TableName[[:blank:]]*=[[:blank:]]*'([^']+)'",REG_EXTENDED);
  r_status&=regcomp(&RE_TableTitle,"TableTitle[[:blank:]]*=[[:blank:]]*'([^']+)'",REG_EXTENDED);
  r_status&=regcomp(&RE_TableType,"TableType[[:blank:]]*=[[:blank:]]*'([^']+)'",REG_EXTENDED);
  r_status&=regcomp(&SBkeys,"![^!][[:alpha:]]",REG_EXTENDED);
  r_status&=regcomp(&SBkey,"(![[:alpha:]][[:alnum:]]*|>([[:alpha:]][[:word:]]*:)*([[:alpha:]][[:word:]])*)",REG_EXTENDED);
  r_status&=regcomp(&SBlink,">([[:alpha:]][[:word:]]*:)*(([[:alpha:]][[:word:]])*)$",REG_EXTENDED);
  
  if (r_status==0){
    printf("regular expressions created.\n");
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
	i_sep=strcspn(s,";,\t");
	fs=s[i_sep]; //field separator
	if (regexec(&RE_TableName,s,2,match,0)==0){
	  TableName=g_strndup(s+match[1].rm_so,match[1].rm_eo-match[1].rm_so);
	  printf("TableName: %s\n",TableName); fflush(stdout);
	}
	if (regexec(&RE_TableType,s,2,match,0)==0){
	  TableType=g_strndup(s+match[1].rm_so,match[1].rm_eo-match[1].rm_so);
	  printf("TableType: %s\n",TableType); fflush(stdout);
	}
	if (regexec(&RE_TableTitle,s,2,match,0)==0){
	  TableTitle=g_strndup(s+match[1].rm_so,match[1].rm_eo-match[1].rm_so);
	  printf("TableTitle: %s\n",TableTitle); fflush(stdout);
	}
      } else if (regexec(&SBkeys,s,1,match,0)==0){
	keys=g_strsplit_set(s,fs,-1);
	int k=g_strv_length(keys);
	for (i=0;i<k;i++) {
	  if (regexec(&SBlink,keys[i],2,match,0)==0){
	    stem=g_strndup(keys[i]+match[1].rm_so,match[1].rm_eo-match[1].rm_so);
	    leaf=g_strndup(keys[i]+match[2].rm_so,match[2].rm_eo-match[2].rm_so);
	    printf("[keys] link to table «%s», ID=«%s» found. I will use ID in hash table.\n",stem,leaf);
	  }
	  g_strstrip(keys[i]);
	}
	sbtab=sbtab_alloc(keys);
      } else if (regexec(&EmptyLine,s,0,NULL,0)){
	printf("skipping empty line.\n");
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
