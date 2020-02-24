#include "event.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* finds time range sub-vectos of a given event. Given a list of
 * event-times and measurement times, it selects all events that
 * happen before measurement time [0] and creates a sub-vector view
 * for those events; this is repeated for all measurement-time
 * points: event_subvector[i] containsd all events that happen right
 * before measurement_time[i] but after measurement_time[i-1].
 */
int /* success/failure flag */
event_sub_list
(gsl_vector *t,/* the experiment's measurement times*/
 event_t *event)/* create sub-vectors for this event*/
{
  size_t i,j;
  assert(t);
  size_t T=t->size;
  size_t K[T];
  double et;
  assert(event);
  size_t N=event->value->size1;
  size_t M=event->value->size2;
  assert(event->time->size == event->value->size1);
  size_t first_TP_after[N]; // first time point of measurement after an event
  /* find out where events belong in the timeline */
  printf("[%s] number of events: %li\n",__func__,N);
  for (i=0;i<N;i++){
    j=0;
    et=gsl_vector_get(event->time,i);
    while (j<T && et>gsl_vector_get(t,j)) j++;
    first_TP_after[i]=(j<T?j:-1);
  }
  event->value_sub=calloc(T,sizeof(gsl_matrix_view));
  event->val_before_t=calloc(T,sizeof(gsl_matrix*));
  event->time_sub=calloc(T,sizeof(gsl_vector_view));
  event->time_before_t=calloc(T,sizeof(gsl_vector*));
  for (j=0;j<T;j++){
    K[j]=0;
    for (i=0;i<N;i++){
      /* count how many events occur right before time point j*/
      if (first_TP_after[i]==j) K[j]++;
    }
    printf("[%s] %li events occur before TimePoint %li (t=%g)\n",__func__,K[j],j,gsl_vector_get(t,j));
    assert(K[j]<=N);
  }
  size_t k=0;
  for (j=0;j<T;j++){
    if(K[j]>0){
      event->value_sub[j]=gsl_matrix_submatrix(event->value,k,0,K[j],M);
      event->val_before_t[j]=&(event->value_sub[j].matrix);
      event->time_sub[j]=gsl_vector_subvector(event->time,k,K[j]);
      event->time_before_t[j]=&(event->time_sub[j].vector);
    }
    k=K[j];
  }
  return EXIT_SUCCESS;
}

/* event targets are specified by name, this function finds the index. If the `target_name` is found in `list_of_names`, this function returns the index, a negative value otherwise*/
int /* the index of `target_name` (negative non failure to find)*/
event_find_target
(char *target_name, /*a \0 terminated string */
 char **list_of_names, /*a list of strings*/
 size_t num) /*the number of items in the `list_of_names`*/
{
  int i;
  for (i=0;i<num;i++){
    if (strcmp(target_name,list_of_names[i])==0){
      return i;
    }
  }
  return -1;
}


/* this adds an event to the event array, an event block can have
 * multiple times at which it happens, but, if an experiment is
 * subject to more than one kind of event (different targets, times,
 * etc.)  then the array is resized and another block is added
*/
event_t* /*returns pointer to newly added event table */
event_add
(event_list *event, /* the event list to modify*/
 gsl_vector *measurement_time, /* for a given experiment*/
 gsl_vector *event_time, /*array of times (table rows)*/
 gsl_vector_int *event_type, /*array of operation types */
 gsl_vector_int *effect, /* whether this affects input or states*/
 char *event_target, /*string of space separated names (operation targets)*/
 gsl_matrix *event_value) /*value subject to operation to be applied to target*/
{
  size_t n=event->num;
  size_t max=event->max_num;
  size_t increment=2;
  if (n==max){
    event->max_num+=increment;
    event->e = realloc(event->e,event->max_num);
  }
  assert(event_time && event_type && event_target && event_value);
  assert(event_type->size==effect->size && event_value->size2 == effect->size);
  event->e[n]=malloc(sizeof(event_t));
  event->e[n]->time=event_time;
  event->e[n]->type=(op_t*) event_type->data;
  event->e[n]->effect=(effect_t*) effect->data;
  event->e[n]->target_names=event_target;
  event->e[n]->value=event_value;
  event_sub_list(measurement_time,event->e[n]);
  event->num++;
  return event->e[n];
}

event_list_t* event_list_alloc(size_t default_size){
  size_t N=default_size>0?default_size:1;
  event_t **e=malloc(sizeof(event_t*)*N);
  assert(e);
  event_list *el=malloc(sizeof(event_list));
  el->num=0;
  el->max_num=N;
  el->e=e;
  return el;
}

void event_list_free(event_list_t *el){
  if (el){
    for (i=0;i<el->num;i++){
      event=el->e[i];
      if(event){
	if(event->value) gsl_matrix_free(event->value);
	if(event->type) free(event->type);
	if(event->effect) free(event->effect);
	if(event->target_names) free(event->target_names);
	if(event->target) gsl_vector_int_free(event->target);
	if(event->time) gsl_vector_free(time);
	if (event->val_before_t) free(val_before_t);
	if (event->value_sub) free(value_sub);
	if (event->time_sub) free(time_sub);
	free(event);
      }
    }
    free(el);
  }
}

/*this function applies the effect of event `e` to either the parameters `p` or the state variables `y` (and the state sensitivity S=dy/dp). The internal model parameters `k` themselves are not subject to change, but `p` contains the input `u`: `p=cat(exp(k),u)`, and `u` can be subject to events*/
void
event_apply
(single_event *e, /* event to apply*/
 gsl_vector *y, /* state to change, possibly*/
 gsl_vector *p, /* parameter vector: [exp(k),u]*/
 gsl_matrix *S) /* state sensitivity.*/
{
  size_t m,M=e->value->size;
  op_t op;
  effect_t effect;
  gsl_vector *yp; /* y or p, depending on effect*/
  double v;
  double V;
  int i;
  for (m=0;m<M;m++){
    /* apply event*/
    i=gsl_vector_get(e->target,m)
    op=e->type[m];
    effect=e->effect[m];
    v=gsl_vector_get(e->value,m);
    switch(effect){
    case event_affects_state:
      assert(i<y->size);
      V=gsl_vector_get(y,i);
      yp=y;
      break;
    case event_affects_input:
      assert(i<p->size);
      V=gsl_vector_get(p,i);
      yp=p;
      break;
    default: yp=NULL;
    }
    if (yp){
      switch(op){
      case event_set: V=v; break;
      case event_add: V+=v; break;
      case event_sub: V-=v; break;
      case event_mul: V*=v; break;
      case event_div: V/=v; break;	
      }
    gsl_vector_set(yp,i,V);
    }
  }
}

/*makes a new single event from a row within an event table*/
event_row_t* /* allocated event */
event_row_link
(size_t point,  /*consider event before this time point*/
 size_t row,    /*row to use*/
 event_t *event)/*table to link event to*/
{
  event_row_t *e=malloc(sizeof(event_row_t));
  e->parent=event;
  e->t=gsl_vector_get(event->time_before_t[point],row);
  e->value_row=gsl_matrix_row(event_table->val_before_t[point],row);
  e->value=&(e->value_row.vector);
  return e;
}


/* this inserts an event table into linked lists the linked lists are
 * timeline global and time ordered; they mix event tables, so
 * different event types can coexist within them.
 * the lists store pointers to the real event table rows.
 */
void events_push
(gsl_vector *time, /* measurement times (data) */
 event_row_t **single, /* array of linked lists */
 event_t *event_table)/* an event table to be inserted into a linked list*/
{
  size_t i,j,J;
  event_row_t *e; /* current event pointer */
  event_row_t *n; /* new event */
  event_row_t **p; /* a pointer to the n's parent.next component */
  
  for (i=0;i<time->size;i++){
    J=event_table->time_before_t[i]->size;
    for (j=0;j<J;j++){
      p=&(single[i]);
      e=*p;
      n=event_row_link(i,j,event_table);
      while (e && n->t < e->t) {
	p=&(e->next);
	e=e->next;
      }
      n->next=e;
      *p=n;
    }
  }
}

/* event_row_links are linked lists*/
size_t /* length of the list*/
get_list_length(event_row_t *s) /* list of event table rows */
{
  size_t l=0;
  while (s) {
    s=s->next;
    l++;
  }
  return l;
}

before_measurement* event_array_alloc(size_t L){
  before_measurement *bt=malloc(sizeof(before_measurement));
  bt->size=L;
  bt->event=malloc(sizeof(event_row_t*)*L);
  return bt;
}

/* this converts a stack of single row event links into a set of arrays for
   quick access. The return value is a list of lists. `before_measurement[i]` refers to all events that happen before measurement time point t[i]. The input `single` is a similar list of stacks: `single[i]` contains event rows prior to t[i]. They are each mapped in reverse order (being stacks) to the return lists.*/
before_measurement** /* a list of arrays `b` with `b[i]` containing
			only events prior to `t[i]` */
convert_to_arrays
(size_t T, /* number of measurement time points */
 event_row_t **single) /* a list of stacks (with length `T`) */
{
  size_t j,L;
  int i;
  before_measurement **b=malloc(sizeof(before_measurement*)*T);
  event_row_t *p;
  if (single){
    for (j=0;j<T;j++){
      L=get_list_length(single[j]);    
      b[j]=event_array_alloc(L);
      p=single[j];
      i=L;
      while (p && i>=0){
	b[j]->event[i--]=p;
	p=p->next
      }
    }
  } else {
    free(b);
    b=NULL;
  }
  return b;
}
