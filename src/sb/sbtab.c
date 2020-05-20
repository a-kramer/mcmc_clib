#include "sbtab.h"
#include <assert.h>
#include "re.h"
#include <math.h>
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
    printf("[%s] received %i keys.\n",__func__,n);
    fflush(stdout);
    for (i=0;i<n;i++) {      
      g_strstrip(keys[i]);
      printf("\t\t\t key(%i) '%s'.\n",i,keys[i]);
      fflush(stdout);
    }    
    sb=malloc(sizeof(sbtab_t));                    // allocate the sbtab_t element
    sb->key=keys;                                  // store the keys as column labels
    sb->column=malloc(sizeof(GPtrArray*)*n);
    //    if (g_strcmp0(keys[0],"!ID")==0 || g_strcmp0(keys[0],"!Name")==0){
    sb->row=g_hash_table_new_full(g_str_hash, g_str_equal, g_free, g_free);
    //}else{
    // sb->row=NULL;
    // }
    sb->col=g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL);
    for (i=0;i<n;i++){
      sb->column[i]=g_ptr_array_new_full(DefaultSize, g_free);
      g_hash_table_insert(sb->col,g_strdup(keys[i]),sb->column[i]);
    }
  }
  printf("[%s] done.\n",__func__);
  fflush(stdout);
  return sb;
}

int sbtab_get_row_index(const sbtab_t* sbtab, const gchar *ID){
  assert(sbtab && sbtab->row);
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
  assert(sbtab_hash);
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
  printf("[%s] Looking up any of: %s.\n",__func__,Names);
  assert(sbtab_hash);
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
  fflush(stdout);
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
      fprintf(stderr,"\n\t\t");
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
  assert(sbtab && sbtab->col);
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
  assert(sbtab && sbtab->col);
  a=g_hash_table_lookup(sbtab->col,key);
  if (a!=NULL) d=g_ptr_array_index(a,i);
  else d=g_strdup("");
  return d;
}

/* given one or more alternative column names, this function looks up the names in the given order until the first match. If no match is found `>` is prepended and the names are tried again.*/
GPtrArray* sbtab_find_column(const sbtab_t *sbtab, const char *key, const char *prefix){
  gchar** Key=g_strsplit(key," ",-1);
  guint n=g_strv_length(Key);
  GPtrArray *a=NULL;
  guint i=0;
  GString *ref_key=g_string_sized_new(strlen(key)+2);
  assert(sbtab && sbtab->col);
  for (i=0;i<n;i++){
    if (prefix) g_string_append(ref_key,prefix);
    g_string_append(ref_key,Key[i]);
    a=g_hash_table_lookup(sbtab->col,ref_key->str);
    if (a) break;
    else ref_key = g_string_assign (ref_key,"");
  }
  /* if (a==NULL){ */
  /*   fprintf(stderr, */
  /* 	    "[%s] warning: lookup of key «%s» in Table «%s» failed.\n", */
  /* 	    __func__,key,sbtab->TableName); */
  /* } */
  g_strfreev(Key);
  g_string_free(ref_key,TRUE);
  return a;
}

int sbtab_find_row(sbtab_t *table, char *id){
  int *row_ptr=NULL;
  int r=-1;
  assert(table && table->row);
  row_ptr=g_hash_table_lookup(table->row,id);
  if (row_ptr) r=row_ptr[0];
  return r;
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


char *MandatoryTableProperty(regex_t *RE,char *s,const char *prop){
  char error_buffer[128];
  regmatch_t match[2];
  char *TableProperty;
  int EC=regexec(RE,s,2,match,0);
  if (EC==0){
    TableProperty=dup_match(&match[1],s);
    printf("Table%s: «%s»\n",prop,TableProperty); fflush(stdout);
  } else {
    regerror(EC,RE,error_buffer,128);
    perror(error_buffer);
    fprintf(stderr,"[%s] error: Table%s is missing.\n",__func__,prop);
    abort();
  }
  return TableProperty;  
}

void print_keys(gchar **keys){
  assert(keys);
  int i,k=g_strv_length(keys);
  printf("[%s] %i headers:",__func__,k); fflush(stdout);
  for (i=0;i<k;i++) {
    assert(keys[i]);
    //g_strstrip(keys[i]);
    printf("«%s» ",keys[i]);
    fflush(stdout);
  }
  printf("done.\n");
  fflush(stdout);
}

sbtab_t* sbtab_from_tsv(char *tsv_file){
  FILE* fid;
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
  int sbtab_header_read=FALSE;
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
	  //printf("[%s] comment found:\n",__func__);
	  //printf("\t\t«%s»\n",s);
	  a=match[0].rm_so;
	  s[a]='\0';
	  g_strstrip(s);
	  //printf("\t\t«%s»\n",s);
	}
	if (sbtab_header_read){
	  if (regexec(&SBkeys,s,2,match,0)==0){
	    keys=g_strsplit_set(s,fs,-1);
	    print_keys(keys);
	    sbtab=sbtab_alloc(keys);
	  } else if (regexec(&EmptyLine,s,0,NULL,0)==0){
	    printf("[%s] skipping empty line «%s».\n",__func__,s);
	  } else {
	    assert(sbtab);
	    sbtab_append_row(sbtab,s,fs);
	  }
	} else {
	  if (regexec(&SBtab,s,0,NULL,0)==0){
	    sbtab_header_read=TRUE;
	    TableName=MandatoryTableProperty(&RE_TableName,s,"Name");
	    TableType=MandatoryTableProperty(&RE_TableType,s,"Type");
	    TableTitle=MandatoryTableProperty(&RE_TableTitle,s,"Title");
	  }
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

/* This is a generic function to convert a column to a numeric matrix
 * Mor specific variants of this exist for some Quantities such as
 * time-vectors and data-matrices. This one allows missing values.
 */
gsl_matrix* sbtab_columns_to_gsl_matrix(sbtab_t *table, GPtrArray *column_names, char *prefix, double default_value){
  int size[2];
  size[0]=table->column[0]->len;
  size[1]=column_names->len;
  gsl_matrix *m=gsl_matrix_alloc(size[0],size[1]);
  gsl_matrix_set_all(m,default_value);
  sbtab_update_gsl_matrix(m,table,column_names,prefix);
  return m;
}

void sbtab_update_gsl_matrix(gsl_matrix *m, sbtab_t *table, GPtrArray *column_names, char *prefix){
  assert(m);
  assert(table);
  assert(column_names);
  int i,j;
  gchar *s,*r;
  gchar *c;
  double value;
  /* figure out the sizes*/
  int size[2];
  size[0]=table->column[0]->len;
  size[1]=column_names->len;
  assert(m->size1 == size[0] && m->size2 == size[1]);
  
  gsl_vector_view m_column;
  GPtrArray *table_column;
  for (i=0;i<size[1]; i++){
    m_column=gsl_matrix_column(m,i);
    c=g_ptr_array_index(column_names,i);
    table_column=sbtab_find_column(table,c,prefix);
    if (table_column){
      for (j=0;j<size[0];j++){
	s=g_ptr_array_index(table_column,j);
	r=s;
	value=strtod(s,&r);
	if (s!=r)  gsl_vector_set(&(m_column.vector),j,value);
      }
    }
  }
}

/* This is a generic function to convert a column to a numeric vector
 * Mor specific variants of this exist for some Quantities such as
 * time-vectors and data-matrices. This one allows missing values.
 */
gsl_vector* sbtab_column_to_gsl_vector(sbtab_t *table,gchar *column_name){
  int j,n;
  GPtrArray *column;
  gsl_vector *v;
  gchar *s,*r;
  double value;
  column=sbtab_find_column(table,column_name,NULL);
  if (column){
    n=column->len;
    v=gsl_vector_alloc(n);
    for (j=0;j<n;j++){
      s=g_ptr_array_index(column,j);
      r=s;
      value=strtod(s,&r);
      if (s==r) value=NAN;
      gsl_vector_set(v,j,value);
    }
  }else{
    v=NULL;
  }
  return v;
}

/* This is a generic function to update a vector (a default)
 * More specific values of this exist for some Quantities such as
 * time-vectors and data-matrices. This one allows missing values.
 */
void sbtab_update_gsl_vector(gsl_vector *v, sbtab_t *table,gchar *column_name){
  int j,n;
  GPtrArray *column;
  gchar *s,*r;
  double value;
  column=sbtab_find_column(table,column_name,NULL);
  if (column){
    n=column->len;
    v=gsl_vector_alloc(n);
    for (j=0;j<n;j++){
      s=g_ptr_array_index(column,j);
      r=s;
      value=strtod(s,&r);
      if (s!=r) gsl_vector_set(v,j,value);
    }
  }
}

