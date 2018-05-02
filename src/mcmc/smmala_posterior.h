int LogPosterior(const double beta, const gsl_vector* x,  void* model_params, double* fx, gsl_vector **dfx, gsl_matrix **FI);
int display_prior_information(const void *Prior_str);
