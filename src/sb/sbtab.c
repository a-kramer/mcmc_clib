#include "sbtab.h"

  SBkey_t *col_label;
  flex_double **column;
  flex_double *row_index;
  GHashTable *row;
  GHashTable *col;


/* bsbtab_t* sbtab_alloc(SBkey_t title) 
 *
 * creates a new table, with column keys and variable length columns;
 * the number of colums is taken from the available keys in the column
 * title list
 */
sbtab_t* sbtab_alloc(SBkey_t *title){
  sbtab_t *sb;
  int i,n;
  char *s;
  assert(title!=NULL);
  n=title->numel;
  sb=malloc(sizeof(sbtab_t));
  sb->col_label=title; // store the keys as column labels
  sb->row=g_hash_table_new_full(g_str_hash, g_str_equal, g_free, flex_array_free);
  sb->col=g_hash_table_new_full(g_str_hash, g_str_equal, g_free, flex_array_free);
  sb->row_index=flex_array_alloc(40);
  sb->column=malloc(sizeof(flex_double*)*n);
  for (i=0;i<n;i++){
    sb->column[i]=flex_array_alloc(40);
    s=sb_key_retrieve(title,i);
    g_hash_table_insert(sb->col,g_strdup(s),sb->column[i]);
  }
  return sb;
}
int sbtab_set(char *key, size_t i){
  return EXIT_SUCCESS;
}
double sbtab_get(char *key, size_t i){
  double d;
  return d;
}
int sbtab_free(){
  return EXIT_SUCCESS;
}
