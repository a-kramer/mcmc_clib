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
#include "sbtab.h"

#define DEFAULT_STR_LENGTH 128

sbtab_t* parse_sb_tab(char *);
 
int main(int argc, char*argv[]){
  int i;
  regex_t *sb_tab_file;
  sbtab_t *table;
  GPtrArray *sbtab;
  char *sptr;
  
  sptr=malloc(sizeof(char)*64);
  sb_tab_file=malloc(sizeof(regex_t));
  regcomp(sb_tab_file,".*[.][tc]sv",REG_EXTENDED);
  sbtab=g_ptr_array_new_full(3,sbtab_free);
  for (i=1;i<argc;i++){
    if (regexec(sb_tab_file,argv[i],0,NULL,0)==0) {
      table=parse_sb_tab(argv[i]);
      g_ptr_array_add(sbtab,table);
    }
    else printf("unknown option «%s».\n",argv[i]);
  }
  guint n=sbtab->len;
  printf("[main] %i tables read.\n",n);
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


sbtab_t* parse_sb_tab(char *sb_tab_file){
  FILE* fid;
  int i;
  char *s;
  char *key;
  size_t n_s=DEFAULT_STR_LENGTH, m_s, length=0;
  regex_t SBtab;
  regex_t RE_TableName, RE_TableTitle, RE_TableType, SBcomment, SBkeys, SBkey;
  gchar *TableName, *TableTitle, *TableType;
  regmatch_t match[2];
  regoff_t a,b;
  int r_status=0;
  gchar **keys;
  sbtab_t *sbtab;
  sbtab=NULL;
  s=malloc(sizeof(char)*n_s);
  r_status&=regcomp(&SBcomment,"[%#]",REG_EXTENDED);
  r_status&=regcomp(&SBtab,"!!SBtab",REG_EXTENDED);  
  r_status&=regcomp(&RE_TableName,"TableName[[:blank:]]*=[[:blank:]]*'([^']+)'",REG_EXTENDED);
  r_status&=regcomp(&RE_TableTitle,"TableTitle[[:blank:]]*=[[:blank:]]*'([^']+)'",REG_EXTENDED);
  r_status&=regcomp(&RE_TableType,"TableType[[:blank:]]*=[[:blank:]]*'([^']+)'",REG_EXTENDED);
  r_status&=regcomp(&SBkeys,"![^!][[:alpha:]]",REG_EXTENDED);
  r_status&=regcomp(&SBkey,"(![[:alpha:]][[:alnum:]]*|>[[:alpha:]][:[:alnum:]]*)",REG_EXTENDED);
  
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
	keys=g_strsplit_set(s,",;\t",-1);
	int k=g_strv_length(keys);
	for (i=0;i<k;i++) g_strstrip(keys[i]);
	sbtab=sbtab_alloc(keys);
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
