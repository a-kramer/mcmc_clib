#include "mcmc_kernel.h"
int write_resume_state(const char *file_name, int rank, int R, const mcmc_kernel *kernel);
int load_resume_state(const char *file_name, int rank, int R, const mcmc_kernel *kernel);
