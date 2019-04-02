#include <stdlib.h>
#include <string.h>

char* flatten(const char **a, size_t n, char *sep){
  int i;
  size_t la=0; // overall length
  size_t nsep=(sep!=NULL?strlen(sep):0);
  size_t l[n];
  for (i=0;i<n;i++){
    l[i]=strlen(a[i]);
    la+=l[i];
  }
  char *s;
  s=malloc(sizeof(char)*(la+nsep*n+1));
  char *to=s;
  for (i=0;i<n;i++){
    to=mempcpy(to,a[i],l[i]);
    if (nsep>0) to=mempcpy(to,sep,nsep);
  }
  to[0]='\0';
  return s;
}
