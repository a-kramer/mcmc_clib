#include "sbtab.h"
#include <assert.h>
#include "re.h"
/* initially, we pre-allocate this:*/
#define DEFAULT_STR_LENGTH 128
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
  //  GPtrArray *C;
  //  C=sbtab_get_column(sbtab,ID);
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
  int r=sbtab_get_row_index(T,D->rowID);
  if (r>=0){
    D->table=T;
  }
}

sbtab_t *sbtab_find_table_with(GHashTable *sbtab_hash, gchar *rowID){
  sbtab_t* table_has_rowID=NULL;
  sbtab_and_row data;
  data.rowID=rowID;
  g_hash_table_foreach(sbtab_hash, sbtab_contains,
		       &data);
  table_has_rowID=data.table;
  if (table_has_rowID){
    printf("[%s] found a table named «%s» that contains «%s».\n",__func__,
	   table_has_rowID->TableName,
	   rowID);
  }else{
    fprintf(stderr,"[%s] warning, no table containing ID «%s» was found.\n",__func__,rowID);
  }
  return table_has_rowID;
}

/* given a list of table names (space separated), this function returns the first table that is found using the given names (in the order given)*/
sbtab_t* sbtab_find(GHashTable *sbtab_hash, const gchar *Names){
  sbtab_t *table=NULL;
  gchar **Name;
  printf("[sbtab_find] Looking up any of: %s.\n",Names);
  Name=g_strsplit(Names," ",-1);
  guint n=(int) g_strv_length(Name);
  //printf("[sbtab_find] %i tokens.\n",n); fflush(stdout);
  gint i=-1;
  for (i=0;i<n;i++){
    table=g_hash_table_lookup(sbtab_hash,Name[i]);
    if (table) break;
  }
  if (table){
    printf("[%s] found table «%s».\n",__func__,table->TableName);
  } else {
    printf("[%s] failed.\n",__func__);
  }
  g_strfreev(Name);  
  return table;
}

/* given a list of space separated table names, this function returns all tables. All names must be valid. Tables that don't exist are not added to the return list. Do not rely on the index in the return list if a table was not found. */
GPtrArray* sbtab_get_tables(GHashTable *sbtab_hash, const gchar *TableNames){
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
      fprintf(stderr,"[%s] warning: %s not found.\n",__func__,Name[i]);
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
      fprintf(stderr,"\n");
      for (i=0;i<n;i++) fprintf(stderr,"«%s» ",s[i]);
      fprintf(stderr,"\n");
      fprintf(stderr,"[%s] input is a tsv file, so delimited by «%s» and inline comments marked by %%\n",__func__,fs);
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

/*returns row `i` from table `sbtab` with column header `key` or empty
  string if not found*/
gchar* sbtab_get(const sbtab_t *sbtab, const char *key, const size_t i){
  gchar *d;
  GPtrArray *a=NULL;
  a=g_hash_table_lookup(sbtab->col,key);
  if (a!=NULL) d=g_ptr_array_index(a,i);
  else d=g_strdup("");
  return d;
}

/* given one or more alternative column names, this function looks up the names in the given order until the first match. If no match is found `>` is prepended and the names are tried again.*/
GPtrArray* sbtab_get_column(const sbtab_t *sbtab, const char *key){
  gchar** Key=g_strsplit(key," ",-1);
  guint n=g_strv_length(Key);
  GPtrArray *a=NULL;
  guint i=0;
  GString *ref_key=g_string_sized_new(strlen(key)+2);
  for (i=0;i<n;i++){
    a=g_hash_table_lookup(sbtab->col,Key[i]);
    if (a) break;
  }
  if(!a){
    fprintf(stderr,
	    "[%s] «%s» not found in «%s» trying keys as references: '>[…]'.\n",
	    __func__,key,sbtab->TableName);
    for (i=0;i<n;i++){
      g_string_printf(ref_key,">%s",Key[i]);
      a=g_hash_table_lookup(sbtab->col,ref_key->str);
      if(a){
	fprintf(stderr,
		"[%s] warning: found «%s» as a reference «%s».\n",
		__func__,Key[i],ref_key->str);
	break;
      }
    }
  }
  if (a==NULL){
    fprintf(stderr,
	    "[%s] warning: lookup of key «%s» (and '>key') in Table «%s» failed.\n",
	    __func__,key,sbtab->TableName);
  }
  g_strfreev(Key);
  g_string_free(ref_key,TRUE);
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


sbtab_t* sbtab_from_tsv(char *tsv_file){
  FILE* fid;
  int i;
  char *s; // string to hold read in lines
  size_t n_s=DEFAULT_STR_LENGTH;
  ssize_t m_s;
  regex_t SBtab;
  regex_t RE_TableName, RE_TableTitle, RE_TableType, SBcomment, SBkeys, SBkey;
  regex_t EmptyLine;
  gchar *TableName, *TableTitle, *TableType;
  regmatch_t match[4];
  regoff_t a;
  int r_status=0;
  gchar **keys;
  sbtab_t *sbtab=NULL;
  
  gchar fs[]="\t";
  s=malloc(sizeof(char)*n_s); // a buffer to read strings from file via getline
  r_status|=egrep(&EmptyLine,"^[[:blank:]]*$");
  r_status|=egrep(&SBcomment,"[%#]");
  r_status|=egrep(&SBtab,"!!SBtab");  
  r_status|=egrep(&RE_TableName,"TableName[[:blank:]]*=[[:blank:]]*'([^']+)'");
  r_status|=egrep(&RE_TableTitle,"TableTitle[[:blank:]]*=[[:blank:]]*'([^']+)'");
  r_status|=egrep(&RE_TableType,"TableType[[:blank:]]*=[[:blank:]]*'([^']+)'");
  r_status|=egrep(&SBkeys,"![^!][[:alpha:]]");
  r_status|=egrep(&SBkey,"(![[:alpha:]][[:alnum:]]*|>([[:alpha:]][[:alnum:]_]*:)*([[:alpha:]][[:alnum:]_])*)");
  assert(r_status==0);
  fid=fopen(tsv_file,"r");
  if (fid==NULL){
    fprintf(stderr,"[%s] file not found «%s».\n",__func__,tsv_file);
    abort();
  }else{
    while (!feof(fid)){
      m_s = getline(&s,&n_s,fid);
      //printf("[%s] %i characters read.\n",__func__,m_s);
      if (m_s>0){
	if (regexec(&SBcomment,s,1,match,0)==0){
	  // remove comment from line
	  a=match[0].rm_so;
	  s[a]='\0';
	}
	if (regexec(&SBtab,s,0,NULL,0)==0){
	  if (regexec(&RE_TableName,s,2,match,0)==0){
	    TableName=dup_match(&match[1],s);
	    printf("TableName: «%s»\n",TableName); fflush(stdout);
	  } else {
	    fprintf(stderr,"[%s] error: TableName is missing.\n",__func__);
	    exit(-1);
	  }
	  if (regexec(&RE_TableType,s,2,match,0)==0){
	    TableType=dup_match(&match[1],s);
	    printf("TableType: «%s»\n",TableType); fflush(stdout);
	  }else {
	    fprintf(stderr,"[%s] warning: TableType is missing.\n",__func__);
	  }
	  if (regexec(&RE_TableTitle,s,2,match,0)==0){
	    TableTitle=dup_match(&match[1],s);
	    printf("TableTitle: «%s»\n",TableTitle); fflush(stdout);
	  }else {
	    fprintf(stderr,"[%s] warning: TableTitle is missing.\n",__func__);
	  }	
	} else if (regexec(&SBkeys,s,2,match,0)==0){
	  keys=g_strsplit_set(s,fs,-1);
	  int k=g_strv_length(keys);
	  //print all headers
	  printf("[%s] %i headers:",__func__,k); fflush(stdout);
	  for (i=0;i<k;i++) {
	    g_strstrip(keys[i]);
	    if (keys[i]) printf("«%s» ",keys[i]);
	    else fprintf(stderr,"key[%i/%i] is NULL. ",i,k); 
	  }
	  printf("done.\n");
	  sbtab=sbtab_alloc(keys);
	} else if (regexec(&EmptyLine,s,0,NULL,0)==0){
	  printf("[%s] skipping empty line «%s».\n",__func__,s);
	} else {
	  assert(sbtab);
	  sbtab_append_row(sbtab,s,fs);
	}
      }
    }
    if (sbtab){
      sbtab->TableTitle=TableTitle;
      sbtab->TableName=TableName;
      sbtab->TableType=TableType;
    }
  }  
  return sbtab;
}
