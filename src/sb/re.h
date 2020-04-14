#ifndef RE_H
#define RE_H
#include <stdlib.h>
#include <stdio.h>
#include <glib.h>
#include <regex.h>

gchar* dup_match(regmatch_t *match, char *source);
int egrep_i(regex_t *RE, char *pattern);
int egrep(regex_t *RE, char *pattern);
void re_match_free(gpointer data);

GPtrArray* ReMatch(regex_t *RE, const char *cs, int num_subexp);

#endif

