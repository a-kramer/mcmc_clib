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
