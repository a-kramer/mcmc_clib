#include <stdlib.h>
#include <stdio.h>
#include <glib.h>

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
  gchar **key;
  //  SBkey_t *col_label;
  // flex_double *col_value;
  // flex_double *row_index;
  GPtrArray **column; //  instead of the hash table col;
  GHashTable *row; // 
  GHashTable *col; //
} sbtab_t;

sbtab_t* sbtab_alloc(gchar **keys);
int sbtab_append_row(const sbtab_t *sbtab, const gchar *data, const char *fs);
GPtrArray* sbtab_get_column(const sbtab_t *sbtab, char *key);
void sbtab_free(void *tab);
sbtab_t* sbtab_find(GHashTable *sbtab_hash, const gchar *Names);
