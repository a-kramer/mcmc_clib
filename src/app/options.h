#include <stdlib.h>
#include <stdio.h>
#ifndef OPTIONS_H
#define OPTIONS_H


typedef struct {
  size_t size;
  size_t max;
  double *value;
} vec_d;

/* collects most of the options that have defaults, values from
   possible configuration files and values from command line
   arguments.
 */
typedef struct {
  char *library_file; /*@code .so@ file*/
  char *output_file; /*name suffix of the result file*/
  char *resume_file; /*name prefix of the sesume file*/
  double target_acceptance; /*taget acceptance of the mcmc algorithm, Metropolis works well with 24%, smmala probably at 50%*/
  double initial_stepsize; /*mcmc step size before tuning*/
  double initial_stepsize_rank_factor; /* the initial size can be
					  adjusted depending on
					  temperature by setting thos
					  to a value larger than 1:
					  @code step_size =
					  pow(rank_factor,rank) *
					  initial_step_size;@*/
  long sample_size; /*target recorded sample size*/
  vec_d *abs_tol; /* ode solver parameter*/
  double rel_tol; /* ode solver parameter*/
  double t0; /* global initial time for integration (ivp)*/
} main_options;  // these are user supplied options to the program

main_options* default_options(char *global_sample_filename_stem,char *lib_name);
void read_abs_tol(vec_d *a, char *abs_tol);


#endif
