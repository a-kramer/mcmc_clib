#ifndef EVENT_GLIB_H
#define EVENT_GLIB_H
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <glib.h>
#include "h5block.h"
#include "../app/event.h"
#include "data.h"
#include <assert.h>
/* this structure holds the properties of an event, which are known once an avent table is opened;*/
typedef struct {
  GArray *target;
  GArray *Op;
  GArray *effect;
  GArray *value;
  GPtrArray *gp_target_name;
  gsl_vector *time;
} g_event_t;

typedef struct {
  GHashTable *sbtab_hash;
  GHashTable *Operation;
  gchar **OP_LABEL;
  GString *EventName;
  gchar *ExperimentID;
  gchar *ExperimentName;
  int ExperimentIndex;
  int ExperimentMajorIndex;
  int ExperimentMinorIndex;
  g_event_t *event;
  data_t *Data;

  //  GPtrArray *event_tab;
  h5block_t *h5;
} EventEnvironment;

g_event_t* g_event_alloc(gsl_vector *t,guint m);
void g_event_free(g_event_t *g_event);
#endif
