#ifndef RE_H
#define RE_H
#include <stdlib.h>
#include <glib.h>
#include <regex.h>

gchar* dup_match(regmatch_t *match, char *source);
int egrep_i(regex_t *RE, char *pattern);
int egrep(regex_t *RE, char *pattern);
#endif
