#ifndef SBTAB_H
#define SBTAB_H
#include <stdlib.h>
#include <stdio.h>
#include <glib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
// "length" means: number of rows, for a table
#define table_length(Table) ((Table && Table->column) ? Table->column[0]->len:0)

/* This struct stores a complete sbtab table.
 * Table{Name,Type,Title} are strings (*gchar)
 *    key: stores a list of column labels (array of strings);
 * column: given a column index, this can be used to retrieve data 
 *    row: find the index of some row given a row ID string, a hash table;
 *         returns pointers to guint values;
 *    col: find the column (ptr-array) given a table column label string (key)
 *         g_hash_table_lookup(col,key[i]) returns a pointer to column[i];
 */
typedef struct {
  gchar *TableTitle;
  gchar *TableName;
  gchar *TableType;
  gchar **key; /* g_strv type of string array, NULL terminated */
  GPtrArray **column; //  instead of the hash table col;
  GHashTable *row; // 
  GHashTable *col; //
} sbtab_t;

sbtab_t* sbtab_from_tsv(char *tsv_file);
sbtab_t* sbtab_alloc(gchar **keys);
int sbtab_append_row(const sbtab_t *sbtab, const gchar *data, const char *fs);
GPtrArray* sbtab_get_column(const sbtab_t *sbtab, const char *key);
void sbtab_free(void *tab);
sbtab_t* sbtab_find(GHashTable *sbtab_hash, const gchar *Names);
GPtrArray* sbtab_get_tables(GHashTable *sbtab_hash, const gchar *TableNames);
int sbtab_get_row_index(const sbtab_t* sbtab, const gchar *ID);
char *sbtab_get_field_by_rowID(const sbtab_t* sbtab, const gchar *ID, const gchar *RequestedColumn);
gchar* sbtab_get(const sbtab_t *sbtab, const char *key, const size_t i);
sbtab_t *sbtab_find_table_with(GHashTable *sbtab_hash, gchar *rowID);
gsl_matrix* sbtab_columns_to_gsl_matrix(sbtab_t *table, GPtrArray *column_names, char *prefix, double default_value);
gsl_vector* sbtab_column_to_gsl_vector(sbtab_t *table, gchar *column_name)
#endif
