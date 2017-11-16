#include <stdlib.h>
#include <regex.h>

typedef struct {
  size_t N;
  size_t numel;
  char **key;
} SBkey_t;

SBkey_t sb_key_alloc(size_t n);
int sb_key_append(SBkey_t *sb, char *token, regmatch_t match);
char * sb_key_retrieve(SBkey_t sb, int i)
