#include "event_glib.h"
/* this structure holds the properties of an event, which are known once an avent table is opened;*/
g_event_t* g_event_alloc(gsl_vector *t,guint m){
  assert(t);
  guint n=t->size;
  g_event_t *g_event=malloc(sizeof(g_event_t));
  g_event->time=t;
  g_event->effect=g_array_sized_new
    (FALSE,
     FALSE,
     sizeof(effect_t), m);
  g_event->Op=g_array_sized_new
    (FALSE,
     FALSE,
     sizeof(op_t), m);
  g_event->target=g_array_sized_new
    (FALSE,
     FALSE,
     sizeof(int), m);    
  g_event->value=g_array_sized_new
    (FALSE,
     FALSE,
     sizeof(double), n*m);
  g_event->gp_target_name=g_ptr_array_sized_new(m);
  return g_event;
}

void g_event_free(g_event_t *g_event){
  g_array_free(g_event->effect,TRUE);
  g_array_free(g_event->Op,TRUE);  
  g_array_free(g_event->target,TRUE);
  g_array_free(g_event->value,TRUE);
  g_ptr_array_free(g_event->gp_target_name,TRUE);
  gsl_vector_free(g_event->time);
}
