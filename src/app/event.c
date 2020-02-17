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
make_sub_events(gsl_vector *t,/* the experiment's measurement times*/
		event_t *event)/* create sub-vectors for this event*/{
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
  event->value_sub=malloc(sizeof(gsl_matrix_view)*T);
  event->val_before_t=malloc(sizeof(gsl_matrix*)*T);
  event->time_sub=malloc(sizeof(gsl_vector_view)*T);
  event->time_before_t=malloc(sizeof(gsl_vector*)*T);
  
  for (j=0;j<T;j++){
    for (i=0;i<N;i++){
      K[j]=0;
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
    }else{
      event->val_before_t[j]=NULL;
      event->time_before_t[j]=NULL;
    }
    k=K[j];
  }
  return EXIT_SUCCESS;
}
 
/* this adds an event to the event array, an event block can have
 * multiple times at which it happens, but, if an experiment is
 * subject to more than one kind of event (different targets, times,
 * etc.)  then the array is resized and another block is added
*/
event_t* /*returns pointer to newly added event table */
add_event(event_table *event,
	  gsl_vector *measurement_time,
	  gsl_vector *event_time, /*array of times (table rows)*/
	  gsl_vector_int *event_type, /*array of operation types (table columns)*/
	  gsl_vector_int *event_target, /*array of operation targets*/
	  gsl_matrix *event_value)/*value subject to
					     operation to be applied
					     to target*/{
  size_t n=event->num;
  size_t max=event->max_num;
  size_t increment=2;
  if (n==max){
    event->max_num+=increment;
    /* the following line re-allocates the event array and returns
       event_t** pointers*/
    event->e = realloc(event->e,event->max_num);
  }
  /* this allocates an actual event and returns an event_t* pointer*/
  event->e[n]=malloc(sizeof(event_t));
  assert(event_time && event_type && event_target && event_value);
  event->e[n]->time=event_time;
  event->e[n]->type=event_type;
  event->e[n]->target=event_target;
  event->e[n]->value=event_value;
  // make measurement time preceeding sub-vectors inside the events
  make_sub_events(measurement_time,event->e[n]);
  event->num++;
  return event->e[n];
}

event_table* event_alloc(size_t default_size){
  size_t N=default_size>0?default_size:1;
  event_t **e=malloc(sizeof(event_t*)*N);
  assert(e);
  event_table *et=malloc(sizeof(event_table));
  et->num=0;
  et->max_num=N;
  et->e=e;
  return et;
}

void event_free(event_t *event){
  if (event){
    gsl_matrix_free(value);
    gsl_vector_int_free(type);
    gsl_vector_int_free(target);
    gsl_vector_free(time);
  }
}

void apply_event(single_event *e, gsl_vector *y, gsl_vector *p){
  size_t m,M=e->value->size;
  int type,target;
  double v;
  double V;
  for (m=0;m<M;m++){
    /* apply event*/
    type=gsl_vector_int_get(e->type,m);
    target=gsl_vector_int_get(e->target,m);
    v=gsl_vector_get(e->value,m);
    if (EVENT_AFFECTS_SPC(type)){
      assert(target<y->size);
      V=gsl_vector_get(y,target);
      switch(type){
      case EVENT_SPC_SET: V=v; break;
      case EVENT_SPC_ADD: V+=v; break;
      case EVENT_SPC_SUB: V-=v; break;
      case EVENT_SPC_MUL: V*=v; break;
      case EVENT_SPC_DIV: V/=v; break;	
      }      
      gsl_vector_set(y,target,V);
    } else if (EVENT_AFFECTS_PAR(type)){
      assert(target<p->size);
      V=gsl_vector_get(p,target);
      switch(type){
      case EVENT_PAR_SET: V=v; break;
      case EVENT_PAR_ADD: V+=v; break;
      case EVENT_PAR_SUB: V-=v; break;
      case EVENT_PAR_MUL: V*=v; break;
      case EVENT_PAR_DIV: V/=v; break;	
      }      
      gsl_vector_set(p,target,V);
    }
  }
}

single_event** single_event_list_alloc(int T){
  single_event **s;
  s=malloc(sizeof(single_event*)*T);
  size_t t;
  for (t=0;t<T;t++) s[t]=NULL;
  return s;
}

/*makes a new single event from a row within an event table*/
single_event* /* allocated event */
new_event(size_t point, /*consider event before this time point*/
	  size_t row, /*row to use*/
	  event_t *event_table)/*table to link event to*/{
  single_event *e=malloc(sizeof(single_event));
  e->t=gsl_vector_get(event_table->time_before_t[point],row);
  e->target=event_table->target;
  e->type=event_table->type;
  e->value_row=gsl_matrix_row(event_table->val_before_t[point],row);
  e->value=&(e->value_row.vector);
  e->next=NULL;
  return e;
}


/* this inserts an event table into linked lists the linked lists are
 * timeline global and time ordered; they mix event tables, so
 * different event types can coexist within them.
 * the lists store pointers to the real event table rows.
 */
void insert_events(gsl_vector *time, /* measurement times (data) */
		   single_event **single, /* array of linked lists */
		   event_t *event_table)/*an event table to be inserted into a linked list*/{
  size_t T=time->size;
  size_t i,j,J;
  single_event *e; /* current event pointer */
  single_event *n; /* new event */
  for (i=0;i<T;i++){
    J=event_table->time_before_t[i]->size;
    e=single[i];
    for (j=0;j<J;j++){
      n=new_event(i,j,event_table);
      if (e){
	if (e->t > n->t){
	  n->next=single[i]; /* insert before root */
	  single[i]=n; /* new root */
	} else {
	  while (e->next!=NULL && (e->next->t < n->t)) e=e->next;
	  /* insert event before next event*/
	  n->next=e->next;
	  e->next=n;	  
	}
      } else { /* create root */
	single[i]=n;
      }
    }
  }
}

size_t get_list_length(single_event *s){
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
  bt->event=malloc(sizeof(single_event*)*L);
  return bt;
}

before_measurement** convert_to_array(size_t T, single_event **single){
  size_t i,j,L;
  before_measurement **before_t=malloc(sizeof(before_measurement*)*T);
  single_event *p;
  if (single){
    for (j=0;j<T;j++){
      L=get_list_length(single[j]);    
      before_t[j]=event_array_alloc(L);
      p=single[j];
      i=0;
      while (p){
	assert(i<L);
	before_t[j]->event[i++]=p;
      }
    }
  } else {
    before_t=NULL;
  }
  return before_t;
}
