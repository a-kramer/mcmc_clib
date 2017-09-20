#define GSL_IS_INT 2
#define GSL_IS_DOUBLE 4
#define GSL_IS_MATRIX 8
#define GSL_IS_VECTOR 16
int gsl_printf(const char *name, void *gsl_thing, int type);
int printf_omp(void *mp);
