#include <stdlib.h>

typedef struct {
  size_t length;
  size_t numel;
  double *data;
} flex_double;

flex_double* flex_array_alloc(size_t n);
int flex_array_set(flex_double *a, double d);
double flex_array_get(flex_double *a, size_t i);
int flex_array_free(flex_double *a);
