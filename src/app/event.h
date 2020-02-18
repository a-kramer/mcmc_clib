#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


/* event types are one of these */
typedef enum effect {event_affects_input, event_affects_state} effect_t;
typedef enum operator {event_set, event_add, event_sub, event_mul, event_div} op_t;
#define EVENT_NIL 0
#define EVENT_SET 1
#define EVENT_ADD 2
#define EVENT_SUB 3
#define EVENT_MUL 4
#define EVENT_DIV 5

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
  gsl_vector_int *type;
  gsl_vector_int *target;
  gsl_matrix *value;
  gsl_matrix **val_before_t; // subvectors of value
  gsl_vector **time_before_t; // sub_vectors of time
  gsl_matrix_view *value_sub; // "views" for the above two 
  gsl_vector_view *time_sub; //
} event_t;

typedef struct {
  size_t num;
  size_t max_num;
  event_t **e;
} event_table;

/* If more than one event table applies to an experiment, the events
 * have to be time ordered. This struct helps to order completely
 * different event types in time.
 */
typedef struct single_event_t single_event;
struct single_event_t {
  double t;
  gsl_vector_int *type;
  gsl_vector_view value_row;
  gsl_vector *value;
  gsl_vector_int *target;
  single_event *next;
};

typedef struct {
  size_t size;
  single_event **event;
} before_measurement;

void insert_events(gsl_vector *time, single_event **single, event_t *event_table);
  
event_t* add_event(event_table *event,
		   gsl_vector *measurement_time,
		   gsl_vector *event_time,
		   gsl_vector_int *event_type, 
		   gsl_vector_int *event_target, 
		   gsl_matrix *event_value);
	
event_table* event_alloc(size_t default_size);
void event_free(event_t *event);
single_event** single_event_list_alloc(int T);
before_measurement** convert_to_array(size_t T, single_event **single);
size_t get_list_length(single_event *s);
void apply_event(single_event *e, gsl_vector *y, gsl_vector *p);
