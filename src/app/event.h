#ifndef EVENT_H
#define EVENT_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/* event types are one of these */
typedef enum effect {event_affects_input, event_affects_state} effect_t;
typedef enum operator {event_set, event_add, event_sub, event_mul, event_div} op_t;
#define NUMBER_OF_EVENTS(event,j) ((event && event->time_before_t[j])?event->time_before_t[j]->size:0)

/* stores event information if the experiment needs that. The plan is:
 * «event» contains a list of sub-vectors with events->time_sub[j]
 * containing all event times i=1:N[j] that happen before
 * t[j]. event->{time,target,type,value} are complete lists of any
 * particular event type. If more than one type of event happens in
 * the experiment, then «event» is resized (it is a pointer array) to
 * accomodate the new event structure.
 */
typedef struct {
  gsl_vector *time;
  op_t *type;             /* operation */
  effect_t *effect;       /* on p or y, effect_t*/
  gsl_vector_int *target; /* index in target vector */
  gsl_matrix *value;
  gsl_matrix **val_before_t;  // subvectors of value
  gsl_vector **time_before_t; // sub_vectors of time
  gsl_matrix_view *value_sub; // "views" for the above two 
  gsl_vector_view *time_sub;  //
} event_t;

/* this is an array of events that influence the same experiment, but
   are different in structure (different number of columns, or
   different column headers)*/
typedef struct {
  size_t num;
  size_t max_num;
  event_t **e;
} event_list_t;

/* If more than one event table applies to an experiment, the events
 * have to be time ordered. This struct helps to order completely
 * different event types in time.
 */
typedef struct event_row event_row_t;
struct event_row {
  double t;
  gsl_vector_view value_row;
  gsl_vector *value;
  event_t *parent;
  event_row_t *next;
};

typedef struct {
  size_t size;
  event_row_t **event;
} before_measurement;

void event_push(event_row_t **single, gsl_vector *time, event_t *event_table);

event_t* event_append(event_list_t *event,
 gsl_vector *measurement_time,
 gsl_vector *event_time,
 gsl_vector_int *event_type,
 gsl_vector_int *effect,
 gsl_vector_int *event_target,
 gsl_matrix *event_value);

void events_are_time_ordered(size_t num, before_measurement **b);
event_list_t* event_list_alloc(size_t default_size);
void event_free(event_t *event);
before_measurement** event_convert_to_arrays(size_t T, event_row_t **single);
size_t list_length(event_row_t *s);
void event_apply(event_row_t *e, gsl_vector *y, gsl_vector *p, gsl_matrix *S);
int event_find_target(char *target_name, const char **list_of_names,  size_t num);
gsl_vector_int* /* the index set of `target_name` (negative non failure to find)*/
event_find_targets
(effect_t *effect, 
 char *target_names, 
 size_t num_targets, 
 const char **list_of_p_names, 
 size_t num_p,
 const char **list_of_x_names,
 size_t num_x);
#endif
