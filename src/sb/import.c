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
gsl_matrix** get_data_matrix(sbtab_t *DataTable,sbtab_t *L1_OUT);

int main(int argc, char*argv[]){
  int i,j;
  int status=EXIT_SUCCESS;
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
      printf("[main] parsing sbtab file %s.\n",argv[i]);
      table=parse_sb_tab(argv[i]);
      g_ptr_array_add(sbtab,table);
      g_hash_table_insert(sbtab_hash,table->TableName,table);
    }
    else printf("unknown option «%s».\n",argv[i]);
  }
  guint n_tab=sbtab->len;
  printf("[main] %i tables read.\n",n_tab);
  // convert all data matrices of numbers into actual double/int numbers
  // 1. find all experiment quantity matrices
  regmatch_t match[5];
  regex_t sb_ListOfExperiments;  
  guint nE;
  guint F;
  gchar *experiment_id,*s;
  gsl_matrix **Y_dY;
  gsl_matrix **Y;
  gsl_matrix **dY;
  double val;
  sbtab_t *DataTable;
  sbtab_t *L1_OUT;
  
  
  // find all output names
  L1_OUT=g_hash_table_lookup(sbtab_hash,"L1_OUT");
  F=L1_OUT->column[0]->len;
  regcomp(&sb_ListOfExperiments,"(([Ll]ist|[Tt]able)Of)?[Ee]xp(eriments)?",REG_EXTENDED);
  
  for (i=0;i<n_tab;i++){
    table=(sbtab_t*) g_ptr_array_index(sbtab,i);
    if (regexec(&sb_ListOfExperiments,table->TableName,5,match,0)==0){
      // table is list of experiments
      nE=table->column[0]->len; // list of experiments -> number of experiments
      printf("[main] found list of experiments: %s with %i enries\n",table->TableName,nE);
      Y=malloc(sizeof(gsl_matrix*)*nE);
      dY=malloc(sizeof(gsl_matrix*)*nE);
      for (j=0;j<nE;j++) {
	experiment_id=(gchar*) g_ptr_array_index(table->column[0],j);
	printf("[main] processing %s\n",experiment_id);
	DataTable=g_hash_table_lookup(sbtab_hash,experiment_id);
	if (DataTable!=NULL){
	  Y_dY=get_data_matrix(DataTable,L1_OUT);
	  Y[j]=Y_dY[0];
	  dY[j]=Y_dY[1];
	} else{
	  fprintf(stderr,"DataTable %s missing (NULL).\n",experiment_id);
	  status&=EXIT_FAILURE;
	}
      }
      printf("result of reading matyrices (Y,dY):\n");
      gsl_matrix_fprintf(stdout,Y[j],"%g,");
    }
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

gsl_matrix** get_data_matrix(sbtab_t *DataTable,sbtab_t *L1_OUT){
  int i,i_r, i_c;
  GPtrArray *c,*dc;
  double val,dval=INFINITY;
  guint F=L1_OUT->column[0]->len;
  guint T=DataTable->column[0]->len;
  gsl_matrix **Y;
  gchar *y,*dy,*s, *ds;
  Y=malloc(sizeof(gsl_matrix*)*2);
  
  for (i=0;i<2;i++) {
    Y[i]=gsl_matrix_alloc(T,F);
  }
  gsl_matrix_set_all(Y[0],1.0);
  gsl_matrix_set_all(Y[1],INFINITY);
  printf("[get_data_matrix] Found experiment %s with %i measurements.\n",DataTable->TableName,T);	  
  for (i_c=0; i_c<F; i_c++){
    y=(gchar *) g_ptr_array_index(L1_OUT->column[0],i_c);
    dy=g_strconcat("SD",y);
    c=sbtab_get_column(y);
    dc=sbtab_get_column(dy);
    if (c!=NULL){
      for (i_r=0; i_r<T; i_r++){
	s = (gchar*) g_ptr_array_index(c,i_r);
	if (dc!=NULL) {
	  ds = (gchar*) g_ptr_array_index(dc,i_r);
	  assert(ds!=NULL);
	  dval=strtod(ds,NULL);
	}
	if (s!=NULL){
	  printf("[get_data_matrix] got «%s» ",s);
	  val=strtod(s,NULL);
	} else {
	  val=1.0;
	  dval=INFINITY;
	}
	//dval=strtod(s,NULL);
	printf(" %g ± %g\n",val,dval);
	gsl_matrix_set(Y[0],i_r,i_c,val);
	gsl_matrix_set(Y[1],i_r,i_c,dval);
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
