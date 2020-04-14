#include "norm.h"
#include "re.h"
#include <assert.h>
#include <string.h>

/* norm_output takes an index `i` and returns an index `j` if output
   `i` is normalised by output `j`. The function's return value is
   `TRUE` if the normalization of output `i` is defined at all,
   `FALSE` otherwise. The value in j[0] remains unchanged if no
   mapping `i->j` is defined. */
int /* TRUE iff mapping `i->j` is defined */
norm_index
(index_map_t M, /* stores normalisation index structure, */
 int i, /* the output index that is perhaps normalised */
 int *out_j) /* stores normalising index `j` upon success */
{
  GArray *N=M;
  int l=N->len;
  assert(0<=i && i<l);
  int j=g_array_index(N,int,i);
  if (j>=0 && j<l){
    assert(out_j);
    out_j[0]=j;
    return TRUE;
  } 
  return FALSE;
}

void norm_index_append(index_map_t M, int j){
  GArray *N=M;
  g_array_append_val(N,j);
}

norm_t* norm_alloc(int nsim, int nout){
  norm_t* N=malloc(sizeof(norm_t));
  assert(nsim>0 && nout>0);
  N->experiment=g_array_sized_new(FALSE,FALSE,sizeof(int),nsim);
  N->time=g_array_sized_new(FALSE,FALSE,sizeof(int),nsim);
  N->output=g_array_sized_new(FALSE,FALSE,sizeof(int),nout);
  return N;
}

typedef struct {
  regex_t *re;
  sbtab_t *table;
  GPtrArray *sub_table;
  GArray *major;
  GArray *minor;
} user_data_t;

user_data_t* user_data_alloc(regex_t *RE, sbtab_t *table, GPtrArray *subt){
  user_data_t *buffer=malloc(sizeof(user_data_t));
  assert(table);
  int n=table_length(table);
  buffer->re=RE;
  buffer->table=table;
  buffer->sub_table=subt;
  buffer->major=g_array_sized_new(FALSE,FALSE,sizeof(int),n);
  if (subt){
    buffer->minor=g_array_sized_new(FALSE,FALSE,sizeof(int),n);
  }
  return buffer;
}

void user_data_free(user_data_t *buffer){
  g_array_free(buffer->major,TRUE);
  if (buffer->sub_table){
    g_array_free(buffer->minor,TRUE);
  }
}

/* g_ptr_array_foreach function*/
void find_reference(gpointer data, gpointer user_data){
  gchar *field=data;
  user_data_t *buffer=user_data;

  GString *Reference;
  int I=UNMAPPED; // any illegal index value will do
  int i=UNMAPPED; // any illegal index value will do
  sbtab_t *st;

  int re_n=(buffer->sub_table==NULL?1:2);
  
  if (field){
    GPtrArray *m=ReMatch(buffer->re,field,re_n);
    if (m->len>=2){
      Reference=g_ptr_array_index(m,1);
      I=sbtab_find_row(buffer->table,Reference->str);
    }
    if (I>=0 && buffer->sub_table && m->len>=3){
      st=g_ptr_array_index(buffer->sub_table,I);
      Reference=g_ptr_array_index(m,2);
      i=sbtab_find_row(st,Reference->str);
    }
    g_ptr_array_free(m,TRUE);
  }
  g_array_append_val(buffer->major,I);
  g_array_append_val(buffer->minor,i);
}

norm_t* normalisation(sbtab_t *ExperimentTable, experiment_type *ExperimentType, sbtab_t *OutputTable, GPtrArray *DataTable, map_t *IdxMap){
  regex_t ID_TimePoint, ID;
  GPtrArray *RelativeTo;

  egrep(&ID_TimePoint,"([[:alpha:]][[:alnum:]_]*)\\[?([^]]+)?\\]?");
  egrep(&ID,"([[:alpha:]][[:alnum:]_]*)");
  
  guint K=IdxMap->major->len;
  int nO=table_length(OutputTable);

  norm_t *N=norm_alloc(K,nO);
  user_data_t *buffer;  

  RelativeTo=sbtab_find_column(OutputTable,"!RelativeTo !NormalisedBy !NormalizedBy",NULL);

  if (RelativeTo){
    buffer=user_data_alloc(&ID,OutputTable,NULL);
    g_ptr_array_foreach(RelativeTo,find_reference, buffer);
    // if this is defined, use g_array_copy
    /* N->output=g_array_copy(buffer->major); */
    // otherwise
    guint l=buffer->major->len;
    N->output=g_array_set_size(N->output,l);
    memcpy(N->output->data,buffer->major->data,sizeof(int)*l);
    
    user_data_free(buffer);
  }
  /* experiment normalisation and time normalisation */
  RelativeTo=sbtab_find_column(ExperimentTable,"!RelativeTo !NormalisedBy !NormalizedBy",NULL);
  int k;
  int I,i;
  int ri,rt;
  int rk,rI;
  if (RelativeTo){
    buffer=user_data_alloc(&ID_TimePoint,ExperimentTable,DataTable);
    g_ptr_array_foreach(RelativeTo,find_reference, buffer);
    for (k=0;k<K;k++){
      I=major_index(IdxMap,k);
      i=minor_index(IdxMap,k);
      ri=i;
      printf("[%s] Experiment %i.%i (%i) normalised by\n",__func__,I,i,k);
      fflush(stdout);
      switch (ExperimentType[I]){
      case dose_response:
	norm_index(buffer->minor,I,&ri);
	rt=0;
	break;
      case time_series:
	ri=0;
	norm_index(buffer->minor,I,&rt);
	break;
      default:
	fprintf(stderr,"[%s] unknown experiment type: %i.\n",__func__,ExperimentType[I]);
	abort();
      }
      rk=UNMAPPED;      
      rI=UNMAPPED;
      if (norm_index(buffer->major,I,&rI) && ISMAPPED(ri)){
	rk=flat_index(IdxMap,rI,ri);
	printf("\texperiment %i.%i (%i flat index)\n",rI,ri,rk);
      } else {
	printf("\tnothing.\n");
      }
      g_array_append_val(N->experiment,rk);
      g_array_append_val(N->time,rt);      
    }
    fflush(stdout);
    user_data_free(buffer);
  }
  regfree(&ID_TimePoint);
  regfree(&ID);
  return N;
}


map_t* empty_map(int nE){
  assert(nE>0);
  map_t *m=malloc(sizeof(map_t));
  int i;
  m->major=g_array_sized_new(FALSE,FALSE,sizeof(int),nE);
  m->minor=g_array_new(FALSE,FALSE,sizeof(int));
  m->flat=g_ptr_array_sized_new(nE);
  for (i=0;i<nE;i++){
    g_ptr_array_add(m->flat,g_array_new(FALSE,FALSE,sizeof(int)));
  }
  return m;
}

void map_free(map_t *m){
  g_array_free(m->major,TRUE);
  g_array_free(m->minor,TRUE);
  GArray *a;
  int i;
  for (i=0;i<m->flat->len;i++){
    a=g_ptr_array_index(m->flat,i);
    g_array_free(a,TRUE);
  }
  g_ptr_array_free(m->flat,FALSE);
  free(m);
}

int major_index(const map_t *m, int flat_index){
  assert(m);
  assert(flat_index < m->major->len);
  int i=g_array_index(m->major,int,flat_index);
  return i;
}

int minor_index(const map_t *m, int flat_index){
  int i;
  assert(m);
  assert(flat_index < m->minor->len);
  i=g_array_index(m->minor,int,flat_index);
  return i;
}

int flat_index(const map_t *m, int major, int minor){
  assert(m && major < m->flat->len);
  const GArray *a=g_ptr_array_index(m->flat,major);
  assert(a && minor < a->len);
  int k=g_array_index(a,int,minor);
  return k;
}


void map_index(map_t *m, int major, int minor){
  assert(major>=0 && minor>=0);
  g_array_append_val(m->major,major);
  g_array_append_val(m->minor,minor);
  GArray *a=g_ptr_array_index(m->flat,major);
  int k=m->major->len;
  g_array_append_val(a,k);
}
