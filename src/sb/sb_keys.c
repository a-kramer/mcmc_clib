#include "sb_keys.h"

/* allocates an SBkey_t object and stores the pointer to a string of
 * all keys in that object. It must already be allocated with malloc
 * SBkey_t sb_key_alloc(char* AllKeys, size_t initial_n);
 *
 * The object of type SBkey_t will then store a list of pointers to 
 * locations where a key starts, after the string has been 
 */
SBkey_t sb_key_alloc(size_t initial_n){
  SBkeys_t *sb;
  sb=malloc(sizeof(SBkeys_t));
  if (sb!=NULL){
    sb->key=malloc(sizeof(char*)*initial_n);
    if (sb->key!=NULL) {
      sb->N=initial_n;
      sb->numel=0;
    }
  }
  return sb;
}

int sb_key_resize(SBkey_t sb, size_t new_n){
  int status=EXIT_SUCCESS;
  if (sb!=NULL && sb->keys!=NULL){
    sb->keys=realloc(sb->keys,sizeof(char*)*new_n);
    if (sb->keys!=NULL) sb->N=new_n;
  } else {
    status &= EXIT_FAILURE;
  }
  return status;
}

/* int sb_key_append(SBkey_t key, match_t match) adds a key from the
 * token string to the list of keys. The key is identified by the
 * offsets stored in the regular expression match_t within the token
 */
int sb_key_append(SBkey_t sb, char *token, match_t *match){
  status=EXIT_SUCCESS;
  size_t m=sb->N+sb->N/10; // potential new size if the key array needs resizing;
  size_t n=sb->numel;      // current number of entries
  regoff_t a=match->rm_so;
  regoff_t b=match->rm_eo;
  regoff_t a,b;
  int l=1+b-a; // length of new key
  m=m>N+1?m:N+1;
  if (n==sb->N) status&=sb_key_resize(sb,m);
  if (n<sb->N) sb->key[n]=calloc(l,sizeof(char));
  if (sb->key[n]!=NULL){
    strncpy(sb->key[n],token+a,l);
    sb->numel++;
  } else status=EXIT_FAILURE;
  return status;
}

char * sb_key_retrieve(SBkey_t sb, int i){
  char *s;
  if (i<sb->numel) s=sb->key[i];
  else s=NULL;
  return s;
}

