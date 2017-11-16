#include <stdlib.h>
#include <stdio.h>
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
#include <glib.h>
#include "flex_array.h"
#include "sb_keys.h"

#define DEFAULT_STR_LENGTH 128

int parse_sb_tab(char *);
 
int main(int argc, char*argv[]){
  int i;
  regex_t *sb_tab_file;
  char *sptr;
  
  sptr=malloc(sizeof(char)*64);
  sb_tab_file=malloc(sizeof(regex_t));
  regcomp(sb_tab_file,".*[.][tc]sv",REG_EXTENDED);

  for (i=1;i<argc;i++){
    if (regexec(sb_tab_file,argv[i],0,NULL,0)==0) parse_sb_tab(argv[i]);
    else printf("unknown option «%s».\n",argv[i]);
  }

  return EXIT_SUCCESS;
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


int parse_sb_tab(char *sb_tab_file){
  FILE* fid;
  char *s;
  char *key;
  size_t n_s=DEFAULT_STR_LENGTH, m_s, length=0;
  regex_t SBtab;
  regex_t TableName, TableTitle, TableType, SBcomment, SBkeys, SBkey;
  char name[DEFAULT_STR_LENGTH];
  regmatch_t match[2];
  regoff_t a,b;
  int r_status=0;
  GHashTable *table;
  
  s=malloc(sizeof(char)*n_s);

  table=g_hash_table_new_full(g_str_hash, g_str_equal, g_free, flex_array_free);
  
  r_status&=regcomp(&SBcomment,"[%#]",REG_EXTENDED);
  r_status&=regcomp(&SBtab,"!!SBtab",REG_EXTENDED);  
  r_status&=regcomp(&TableName,"TableName[[:blank:]]*=[[:blank:]]*'([^']+)'",REG_EXTENDED);
  r_status&=regcomp(&TableTitle,"TableTitle[[:blank:]]*=[[:blank:]]*'([^']+)'",REG_EXTENDED);
  r_status&=regcomp(&TableType,"TableType[[:blank:]]*=[[:blank:]]*'([^']+)'",REG_EXTENDED);
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
	if (regexec(&TableName,s,2,match,0)==0){
	  copy_match(name,&match[1],s,DEFAULT_STR_LENGTH);
	  printf("TableName: %s\n",name); fflush(stdout);
	}
	if (regexec(&TableType,s,2,match,0)==0){
	  copy_match(name,&match[1],s,DEFAULT_STR_LENGTH);
	  printf("TableType: %s\n",name); fflush(stdout);
	}
	if (regexec(&TableTitle,s,2,match,0)==0){
	  copy_match(name,&match[1],s,DEFAULT_STR_LENGTH);
	  printf("TableTitle: %s\n",name); fflush(stdout);
	}
      } else if (regexec(&SBkeys,s,1,match,0)==0){
	SBkey_t *sb;
	int k=0;
	sb=sb_key_alloc(12);
	assert(sb!=NULL);	

	flex_double *column;
	column=flex_array_alloc(60);
	key=strtok(s,",;\t");
	while (key!=NULL){
	  printf("strtok: «%s»,\t",key);
	  if (regexec(&SBkey,key,2,match,0)==0){
	    sb_key_append(sb,key,&match[1]);	    
	    //	    copy_match(name,&match[1],key,DEFAULT_STR_LENGTH);
	    printf("SBkey: «%s»\n",sb_key_retrieve(sb,k));
	    k++;
	    g_hash_table_insert(table,g_strdup(name),);
	  }
	  key=strtok(NULL,",;\t");
	}
      }
      
    }
  }
  g_hash_table_destroy(table);
  return EXIT_SUCCESS;
}
