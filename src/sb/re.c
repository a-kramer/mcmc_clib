#include "re.h"
#include <stdio.h>
#include <errno.h>
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

int egrep_i(regex_t *RE, char *pattern){
  char error_buffer[128];
  int EC=regcomp(RE,pattern,REG_EXTENDED|REG_ICASE);
  if (EC) {
    fprintf(stderr,"[%s] regcomp failed with pattern «%s».\n",__func__,pattern);
    regerror(EC,RE,error_buffer,128);
    perror(error_buffer);
    abort();
  }  
  return EC;
}

int egrep(regex_t *RE, char *pattern){
  char error_buffer[128];
  int EC=regcomp(RE,pattern,REG_EXTENDED);
  if (EC) {
    fprintf(stderr,"[%s] regcomp failed with pattern «%s».\n",__func__,pattern);
    regerror(EC,RE,error_buffer,128);
    perror(error_buffer);
    abort();
  }  
  return EC;
}

void re_match_free(gpointer data){
  g_string_free((GString*) data,TRUE);
}

/*performs a regular expression match and stores all matches in an array of `GString`s */
GPtrArray* /*array of GStrings with matches: 0 is the whole match, i>0 matched subexpressions (all are copies) */
ReMatch(regex_t *RE, /* regular expression to match (compiled) */
 const char *cs, /* string to perform match on */
 int num_subexp) /* number of parethesised subexpressions */
{
  int i;
  int num=num_subexp>0?num_subexp+1:1;
  regmatch_t match[num];
  char error_buffer[128];
  GString *temp;
  GPtrArray *m=g_ptr_array_new_full(num, re_match_free);
  regoff_t a,b;
  int length;

  int EC;
  EC=regexec(RE, cs, num, match, REG_EXTENDED);
  switch (EC){
  case 0:
    for (i=0;i<num;i++){
      temp=g_string_sized_new(32);
      a=match[i].rm_so;
      b=match[i].rm_eo;
      length=b-a;
      g_string_append_len(temp,&(cs[a]),length);
      g_ptr_array_add(m,temp);
    }
    break;
  case REG_NOMATCH:
    /* return empty match array */
    for(i=0;i<num;i++) g_ptr_array_add(m,g_string_sized_new(0));
    fprintf(stderr,"[%s] match failed.\n",__func__);
    break;
  default: 
    regerror(EC,RE,error_buffer,128);
    perror(error_buffer);
    abort();   
  }
  return m;
}
