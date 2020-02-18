#include "sbtab.h"
#include <assert.h>
/* sbtab_t* sbtab_alloc(gchar **keys) 
 *
 * creates a new table given column keys. The columns are stored as
 * variable length arrays (GPtrArray); the number of colums is taken from
 * the list of available keys.
 */
sbtab_t* sbtab_alloc(gchar **keys){
  sbtab_t *sb=NULL;
  int i,n;
  guint DefaultSize=32;
  if (keys!=NULL){
    n=g_strv_length(keys);
    printf("[sbtab_alloc] received %i keys.\n",n);
    for (i=0;i<n;i++) {      
      g_strstrip(keys[i]);
      printf("\t\t\t key(%i) '%s'.\n",i,keys[i]);
    }    
    sb=malloc(sizeof(sbtab_t));                    // allocate the sbtab_t element
    sb->key=keys;                                  // store the keys as column labels
    sb->column=malloc(sizeof(GPtrArray*)*n);
    if (g_strcmp0(keys[0],"!ID")==0 || g_strcmp0(keys[0],"!Name")==0){
      sb->row=g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
    }else{
      sb->row=NULL;
    }
    sb->col=g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
    for (i=0;i<n;i++){
      sb->column[i]=g_ptr_array_new_full(DefaultSize, g_free);
      g_hash_table_insert(sb->col,g_strdup(keys[i]),sb->column[i]);
    }
  }
  return sb;
}

int sbtab_get_row_index(const sbtab_t* sbtab, const gchar *ID){
  assert(sbtab);
  GHashTable *row=sbtab->row;
  guint *r;
  int i=-1;
  r=g_hash_table_lookup(row,ID);
  if (r){
    i=r[0];
  }
  return i;
}

char* sbtab_get_field_by_rowID(const sbtab_t* sbtab, const gchar *ID, const gchar *RequestedColumn){
  assert(sbtab);
  assert(ID);
  int r=sbtab_get_row_index(sbtab,ID);
  GPtrArray *C;
  C=sbtab_get_column(sbtab,ID);
  char *c=sbtab_get(sbtab,ID,r);
  return c;
}

typedef struct {
  sbtab_t *table;
  gchar *rowID;
} sbtab_and_row;

void sbtab_contains(gpointer key, gpointer value, gpointer data){
  sbtab_t *T=value;
  sbtab_and_row *D=data;
  int r=sbtab_get_row_index(T,data->rowID);
  if (r>=0){
    D->table=T;
  }
}

sbtab_t *sbtab_find_table_with(const GHashTable *sbtab_hash, const gchar *rowID){
  sbtab_t* table_has_rowID=NULL;
  sbtab_and_row data;
  data->rowID=rowID;
  g_hash_table_foreach(sbtab_hash, sbtab_contains,
		       &data);
  table_has_rowID=data.table;
  if (table_has_rowID){
    printf("[%s] found a table named «%s» that contains «%s».\n",
	   table_has_rowID->TableName,
	   rowID);
  }else{
    fprintf(stderr,"[%s] warning, no table containing ID «%s» was found.\n",rowID);
  }
  return table_has_rowID;
}

sbtab_t* sbtab_find(const GHashTable *sbtab_hash, const gchar *Names){
  sbtab_t *table=NULL;
  gchar **Name;
  printf("[sbtab_find] Looking up any of: %s.\n",Names);
  Name=g_strsplit(Names," ",-1);
  guint n=(int) g_strv_length(Name);
  //printf("[sbtab_find] %i tokens.\n",n); fflush(stdout);
  gint i=-1;
  while (table==NULL && i+1<n){
    table=g_hash_table_lookup(sbtab_hash,Name[++i]);
  }
  if (table!=NULL){
    printf("[%s] found table «%s».\n",__func__,Name[i]);
  } else {
    printf("[%s] not found (%s).\n",__func__,Name[i]);
  }
  g_strfreev(Name);  
  return table;
}

GPtrArry* sbtab_get_tables(const GHashTable *sbtab_hash, const gchar *TableNames){
  assert(TableNames);
  gchar **Name=g_strsplit(TableNames," ",-1);
  guint n=g_strv_length(Name);
  GPtrArray *table=g_ptr_array_sized_new(n);
  sbtab_t *t;
  guint i;
  for (i=0;i<n;i++){
    t=g_hash_table_lookup(sbtab_hash,Name[i]);
    if (t){
      g_ptr_array_add(table,t);
      printf("[%s] %s found.\n",__func__,Name[i]);
    } else {
      printf("[%s] warning: %s not found.\n",__func__,Name[i]);
    }
  }
  return table;
}

int sbtab_append_row(const sbtab_t *sbtab, const char *data, const char *fs){
  int i,N;
  gchar **s;
  guint *r;
  int status=EXIT_SUCCESS;
  s=g_strsplit_set(data,fs,-1);
  if (data && s){
    int n=g_strv_length(s);
    guint c=g_hash_table_size(sbtab->col);
    N=c<n?c:n;
    for (i=0;i<n;i++) g_strstrip(s[i]);
    assert(sbtab);
    if (sbtab->row){
      // get number of rows;
      r=malloc(sizeof(guint));
      r[0]=g_hash_table_size(sbtab->row);
      g_hash_table_insert(sbtab->row,g_strdup(s[0]),r);
      printf("+ %i «%s»\n",r[0],s[0]);
    }
    for (i=0;i<N;i++){
      g_ptr_array_add(sbtab->column[i],g_strdup(s[i]));
    }      
    if(n!=c) {
      fprintf(stderr,"[%s] data contained %i entries while table has %i headers\n\t\t",__func__,n,c);
      for (i=0;i<c;i++) fprintf(stderr,"«%s» ",sbtab->key[i]);
      fprintf(stderr,"\n\t\t");
      for (i=0;i<n;i++) fprintf(stderr,"«%s» ",s[i]);
      fprintf(stderr,"\t\t\tinput is a [tc]sv file, so delimited by ,; or \\t and inline comments marked by %%\n");
      guint nTK;
      const gchar **TK=(const gchar **) g_hash_table_get_keys_as_array(sbtab->col,&nTK);
      for (i=0;i<c;i++)
	fprintf(stderr,"\t\t\tcolumn(%i): '%s'.\n",i,TK[i]);
      status|=EXIT_FAILURE;
    }
  }
  return status;
}

int sbtab_append(const sbtab_t *sbtab, const char *key, char *data){
  GPtrArray *a;  
  a=g_hash_table_lookup(sbtab->col,key);
  if (a!=NULL) {
    g_ptr_array_add(a,g_strdup(data));
    return EXIT_SUCCESS;
  }else{
    return EXIT_FAILURE;
  }
}

gchar* sbtab_get(const sbtab_t *sbtab, char *key, size_t i){
  gchar *d;
  GPtrArray *a=NULL;
  a=g_hash_table_lookup(sbtab->col,key);
  if (a!=NULL) d=g_ptr_array_index(a,i);
  else d=NULL;
  return d;
}

GPtrArray* sbtab_get_column(const sbtab_t *sbtab, const char *key){
  gchar** Key=g_strsplit(key," ",-1);
  guint n=g_strv_length(Key);
  GPtrArray *a=NULL;
  guint i=0;
  while (!a && i<n){
    a=g_hash_table_lookup(sbtab->col,Key[i++]);
  }
  if (a==NULL){
    fprintf(stderr,"[sbtab_get_column] warning: lookup of «%s» in Table «%s» failed.\n",key,sbtab->TableName);
  }
  g_strv_free(Key);
  return a;
}

void sbtab_free(void *tab){
  sbtab_t *sbtab=tab;
  int n=g_strv_length(sbtab->key);
  int i;
  g_hash_table_destroy(sbtab->col);
  g_hash_table_destroy(sbtab->row);
  for (i=0;i<n;i++){    
    g_ptr_array_free(sbtab->column[i],TRUE);    
  }
  g_free(sbtab->column);
  g_strfreev(sbtab->key);
  free(sbtab);
}
