/*
 * SMMALA implementation
 */
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include <string.h>
#include <mpi.h>
#include "smmala.h"
#include "mv_norm.h"
#include "model_parameters_smmala.h"

#define SWAP(a, b, tmp)  tmp = a; a = b; b = tmp;

// fx = beta*lx + px
typedef struct{
  double fx;			/* store current posterior, including beta */
  double stepsize;	/* store current step size */
  double* dfx;		/* store current likelihood gradient */
  double* Hfx;		/* store current cholesky factor of FI */
  
  double* new_x;		/* tmp working space for new x */
  double* new_dfx;
  double* new_Hfx;
  double* mean_vec;
  double* cholH_mat;
  double target_acceptance;
} smmala_params;

smmala_comm_buffer* smmala_comm_buffer_alloc(int D){
  smmala_comm_buffer *buffer;
  int overall_error=-1;
  
  buffer=malloc(sizeof(smmala_comm_buffer));
  if (buffer==NULL) {
      printf("buffer is a NULL pointer\n");
      MPI_Abort(MPI_COMM_WORLD, overall_error);
  }
  buffer->x=malloc(sizeof(double)*D);
    if (buffer->x==NULL) {
      printf("buffer x is a NULL pointer\n");
      MPI_Abort(MPI_COMM_WORLD, overall_error);
  }
  buffer->lx=malloc(sizeof(double)*1);
  if (buffer->lx==NULL) {
      printf("buffer lx is a NULL pointer\n");
      MPI_Abort(MPI_COMM_WORLD, overall_error);
  }

  buffer->dlx=malloc(sizeof(double)*D);
  if (buffer->dlx==NULL) {
      printf("buffer dlx is a NULL pointer\n");
      MPI_Abort(MPI_COMM_WORLD, overall_error);
  }
  buffer->px=malloc(sizeof(double)*1);
  if (buffer->px==NULL) {
      printf("buffer px is a NULL pointer\n");
      MPI_Abort(MPI_COMM_WORLD, overall_error);
  }
  buffer->dpx=malloc(sizeof(double)*D);
  if (buffer->dpx==NULL) {
      printf("buffer dpx is a NULL pointer\n");
      MPI_Abort(MPI_COMM_WORLD, overall_error);
  }
  buffer->FI_l=malloc(sizeof(double)*D*D);
  if (buffer->FI_l==NULL) {
      printf("buffer FI_l is a NULL pointer\n");
      MPI_Abort(MPI_COMM_WORLD, overall_error);
  }
  buffer->FI_p=malloc(sizeof(double)*D*D);
  if (buffer->FI_p==NULL) {
      printf("buffer FI_p is a NULL pointer\n");
      MPI_Abort(MPI_COMM_WORLD, overall_error);
  }
  return buffer;
}

int smmala_exchange_information(mcmc_kernel* kernel, const int DEST, void *buffer){
  MPI_Status status; // MPI_Status contains: MPI_SOURCE, MPI_TAG, MPI_ERROR
  int N=kernel->N; // size of x, dfx; N*N is size of FI;
  smmala_params* state = (smmala_params*) kernel->kernel_params;
  smmala_model* model = (smmala_model*) kernel->model_function;
  ode_model_parameters *omp = (ode_model_parameters*) model->m_params;  
  int TAG=0;
  int SRC=DEST; // send to and receive from the same process
  smmala_comm_buffer *buf=buffer;
  //error handling for nerr communication steps
  int nerr=7;
  int i, ec=MPI_SUCCESS, overall_error=MPI_SUCCESS;
  int len, error_class;
  char error_string[MPI_MAX_ERROR_STRING];

  //things to be communicated:
  double *send_buffer[]={kernel->x,omp->FI_l->data, omp->prior_inverse_cov->data, &(omp->lx), omp->dlx->data, &(omp->px), omp->dpx->data};
  double *recv_buffer[]={   buf->x,      buf->FI_l,                    buf->FI_p,    buf->lx,       buf->dlx,    buf->px, buf->dpx};
  int            size[]={        N,            N*N,                          N*N,          1,              N,          1,        N};
  //int r;
  //MPI_Comm_rank(MPI_COMM_WORLD,&r);
  //printf("[rank %i] exchange information\n",r);  
  // every kernel needs to know x and fx;
  for (i=0;i<nerr;i++){
    /* if (send_buffer[i]==NULL) { */
    /*   printf("[rank %i] send_buffer[%i] is a NULL pointer\n",r,i); */
    /*   MPI_Abort(MPI_COMM_WORLD, overall_error); */
    /* } */
    /* if (recv_buffer[i]==NULL) { */
    /*   printf("[rank %i] recv_buffer[%i] is a NULL pointer\n",r,i); */
    /*   MPI_Abort(MPI_COMM_WORLD, overall_error); */
    /* } */

    ec=MPI_Sendrecv(   send_buffer[i], size[i], MPI_DOUBLE, DEST, TAG,
		       recv_buffer[i], size[i], MPI_DOUBLE, SRC, TAG,
		       MPI_COMM_WORLD, &status);
    overall_error&=ec;
    overall_error&=status.MPI_ERROR;
    if (ec != MPI_SUCCESS){
      MPI_Error_string(ec, error_string, &len);
      fprintf(stderr, "[Comm with rank %3i, part %i] error: %s\n", DEST, i, error_string);
      
      MPI_Error_class(ec, &error_class);
      MPI_Error_string(error_class, error_string, &len);
      fprintf(stderr, "[Comm with rank %3i, part %i] class: %s\n", DEST, i, error_string);
    }
  }
  //printf("[rank %i] exchange information done\n",r);  
  if (overall_error != MPI_SUCCESS){
   fprintf(stderr, "[Comm with rank %3i] Some Communication error occured.\n", DEST);
   MPI_Error_string(overall_error, error_string, &len);
   fprintf(stderr, "[Comm with rank %3i] %s\n", DEST, error_string);
   MPI_Abort(MPI_COMM_WORLD, overall_error);
  }
  /* FILE *f; */
  /* char name[128]; */
  /* int r; */
  /* MPI_Comm_rank(MPI_COMM_WORLD,&r); */
  /* sprintf(name,"rank_%i_comm_results.txt",r); */
  /* f=fopen(name,"a"); */
  /* MPI_Barrier(MPI_COMM_WORLD); */
  /* fprintf(f,"[rank %i] my x: ",r);   for (i=0;i<N;i++)   fprintf(f," %+g ",kernel->x[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] my fx: ",r);                      fprintf(f," %+g ",kernel->fx[0]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] my dfx: ",r); for (i=0;i<N;i++)   fprintf(f," %+g ",state->dfx[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] my Hfx: ",r); for (i=0;i<N*N;i++) fprintf(f," %+g ",state->Hfx[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] my lx: ",r);                      fprintf(f," %+g ",omp->lx); fprintf(f,"\n"); */
  /* MPI_Barrier(MPI_COMM_WORLD); */
  /* fprintf(f,"[rank %i] recv x: ",r);   for (i=0;i<N;i++)   fprintf(f," %+g ",buf->x[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] recv fx: ",r);                      fprintf(f," %+g ",buf->fx[0]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] recv dfx: ",r); for (i=0;i<N;i++)   fprintf(f," %+g ",buf->dfx[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] recv Hfx: ",r); for (i=0;i<N*N;i++) fprintf(f," %+g ",buf->FI[i]); fprintf(f,"\n"); */
  /* fprintf(f,"[rank %i] recv lx: ",r);                      fprintf(f," %+g ",buf->lx[0]); fprintf(f,"\n"); */
  /* fclose(f); */
  return overall_error;
}

int smmala_swap_chains(mcmc_kernel* kernel, const int master, const int rank, const int other_rank, void *buffer){
  double a;
  int swap_accepted=0;
  int N=kernel->N; // size of x, dfx; N*N is size of FI;
  gsl_rng* rng = (gsl_rng*) kernel->rng;
  smmala_params* state = (smmala_params*) kernel->kernel_params;
  smmala_model* model = (smmala_model*) kernel->model_function;
  ode_model_parameters *omp = (ode_model_parameters*) model->m_params;
  int i;
  int TAG=0;
  smmala_comm_buffer *buf=buffer;
  //  log likelihood of x: omp->lx;
  //  log      prior of x: omp->px;
  double beta=omp->beta;
  double their_beta;
  MPI_Sendrecv(&beta, 1, MPI_DOUBLE, other_rank, TAG, &their_beta, 1, MPI_DOUBLE, other_rank, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  TAG++;
  a=gsl_sf_exp((their_beta - beta)*(omp->lx - buf->lx[0]));
  //printf("[rank %i] swap probability: exp((%+g - %+g)*(%+g - %+g))=%+g\n",rank,their_beta,beta,omp->lx,buf->lx[0],a);
  if (master){
    double r1=gsl_rng_uniform(rng);
    swap_accepted=r1<a?1:0;
    MPI_Send(&swap_accepted, 1, MPI_INT, other_rank,TAG,MPI_COMM_WORLD);
  }else{
    MPI_Recv(&swap_accepted, 1, MPI_INT, other_rank,TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (swap_accepted!=0){ // both sides copy their received comm buffer into the kernel
    memcpy(kernel->x,buf->x,N*sizeof(double));
    omp->lx=buf->lx[0];
    omp->px=buf->px[0];
    memcpy(omp->dpx->data,buf->dpx,sizeof(double)*N);
    memcpy(omp->dlx->data,buf->dlx,sizeof(double)*N);
    memcpy(omp->prior_inverse_cov->data,buf->FI_p,sizeof(double)*N*N);
    memcpy(omp->FI_l->data,buf->FI_l,sizeof(double)*N*N);
    
    kernel->fx[0]=beta*omp->lx + omp->px;

    for (i=0;i<N;i++){
      state->dfx[i]=beta*gsl_vector_get(omp->dlx,i)+gsl_vector_get(omp->dpx,i);
    }
    for (i=0;i<N*N;i++){
      state->Hfx[i]=gsl_pow_2(beta)*omp->FI_l->data[i]+omp->prior_inverse_cov->data[i];
    }    
  }
  return swap_accepted;
}

static smmala_params* smmala_params_alloc(int N, double step_size, double target_acceptance){
	
 smmala_params* params = (smmala_params*) malloc( sizeof(smmala_params) );
  if (params == 0){
    /* TODO: write a proper error handler */
    fprintf(stderr,"smmala_params_alloc: malloc failed to allocate memory for smmala_params \n");
    return 0;
  }

  params->target_acceptance=target_acceptance;

  params->dfx     = 0; params->Hfx		= 0;
  params->new_x   = 0; params->new_dfx	= 0;
  params->new_Hfx = 0; params->mean_vec	= 0;
  params->cholH_mat = 0;
	
  params->dfx = (double*) malloc( N * sizeof(double) );
  if (params->dfx == 0)
	  goto error_cleanup;

  params->Hfx = (double*) malloc( N * N * sizeof(double) );
  if (params->Hfx == 0)
	  goto error_cleanup;

  params->new_x = (double*) malloc( N * sizeof(double) );
  if (params->new_x == 0 ) 
	  goto error_cleanup;
  
  params->new_dfx = (double*) malloc( N * sizeof(double) );
  if (params->new_dfx == 0 ) 
	goto error_cleanup;
	
  params->new_Hfx = (double*) malloc( N * N * sizeof(double) );
  if (params->new_Hfx == 0 ) 
	  goto error_cleanup;

  params->mean_vec = (double*) malloc(  N * sizeof(double) );
  if (params->mean_vec == 0 ) 
	  goto error_cleanup;
	
  params->cholH_mat = (double*) malloc( N * N * sizeof(double) );
  if (params->cholH_mat == 0 ) 
	  goto error_cleanup;
	
  params->fx = 0;
  params->stepsize = step_size;

  return(params);
	
error_cleanup:
	fprintf(stderr," smmala_params_alloc: malloc failed \n");
	
	if (params->dfx)
		free(params->dfx);
	
	if (params->Hfx)
		free(params->Hfx);
	
	if (params->new_x)
		free(params->new_x);
	
	if (params->new_dfx)
		free(params->new_dfx);
	
	if (params->new_Hfx)
		free(params->new_Hfx);

	if (params->mean_vec)
		free(params->mean_vec);
	
	if (params->cholH_mat)
		free(params->cholH_mat);

	free(params);
	return 0;
}

static void smmala_params_free(smmala_params* params){
	free(params->cholH_mat);
	free(params->mean_vec);
	free(params->new_Hfx);
	free(params->new_dfx);
	free(params->new_x);
	free(params->Hfx);
	free(params->dfx);
	free(params);
}



static int smmala_kernel_init(mcmc_kernel* kernel, const double* x){
  int res,i,n;
  n = kernel->N;
	
  smmala_params* params = (smmala_params*) kernel->kernel_params;
  /* copy x to the kernel x state */
  for ( i=0; i < n; i++)
	  kernel->x[i] = x[i];

  smmala_model* model = kernel->model_function;
	
  res = model->Likelihood(x, model->m_params, &(params->fx),
						  params->dfx, params->Hfx );
  if (res != 0){
    fprintf(stderr,"smmala_kernel_init: Likelihood function failed\n");
    exit(-1);
  }
  kernel->fx = &(params->fx);	
  gsl_matrix_view Hfx = gsl_matrix_view_array(params->Hfx,n,n); 
  res = gsl_linalg_cholesky_decomp( &Hfx.matrix );
  if (res != 0){
	fprintf(stderr,"Error: matrix not positive definite in smmala_init.\n");
	return -1;
  }
  
  return 0;
}

static int smmala_kernel_init_rand(mcmc_kernel* kernel){
	int res,i,n;
	n = kernel->N;
	
	gsl_rng* rng = (gsl_rng*) kernel->rng;
	smmala_model* model = kernel->model_function;
	
	/* sample random x from the prior */
	double x[n];		/* automatic alloc */
	model->Prior_rnd(rng, model->m_params, x);
	
	/* copy x to the kernel x state */
	for ( i=0; i < n; i++)
		kernel->x[i] = x[i];
	
	smmala_params* params = (smmala_params*)kernel->kernel_params;
	
	res = model->Likelihood(x, model->m_params, &(params->fx),
							params->dfx, params->Hfx );
	/* TODO: write a proper error handler */
	if (res != 0){
		fprintf(stderr,"smmala_kernel_init: Likelihood function failed\n");
		return 1;
	}
	
	gsl_matrix_view Hfx = gsl_matrix_view_array(params->Hfx,n,n); 
	res = gsl_linalg_cholesky_decomp( &Hfx.matrix );
	if (res != 0){
		fprintf(stderr,"Error: matrix not positive definite in smmala_init.\n");
		return -1;
	}
	
	return 0;
}


static int nat_grad_step(const gsl_vector* x, const gsl_matrix* cholPr, 
						 const gsl_vector* grad, gsl_vector* result, double stepSize){

	/* TODO: error checking and return error code */
	
	/* newx = x + (0.5*e)*H^{-1}*grad */
	/* H^{-1}*grad is similar to solving the system H*a = grad */
	gsl_linalg_cholesky_solve(cholPr, grad, result);
	gsl_vector_scale(result, 0.5*stepSize);
	gsl_vector_add(result, x);
	
	return 0;
}


static int smmala_kernel_sample(mcmc_kernel* kernel, int* acc){
  /* TODO: check return values for errors from gsl functions */

  smmala_params* state = (smmala_params*) kernel->kernel_params;
  smmala_model* model = (smmala_model*) kernel->model_function;
  gsl_rng* rng = (gsl_rng*) kernel->rng;
  
  int n = kernel->N;
  double stepsize = state->stepsize;
  
  /* views are allocated in the stack and will be cleared when they get out of scope */
  gsl_vector_view x_vec_v = gsl_vector_view_array(kernel->x, n);
  gsl_vector_view dfx_vec_v = gsl_vector_view_array(state->dfx, n);
  gsl_matrix_view H_mat_v = gsl_matrix_view_array(state->Hfx, n, n);
  gsl_vector_view mean_vec_v = gsl_vector_view_array(state->mean_vec, n);
  gsl_matrix_view cholH_mat_v = gsl_matrix_view_array(state->cholH_mat, n, n);
  gsl_vector_view new_x_v = gsl_vector_view_array(state->new_x, n);

  gsl_matrix_memcpy(&cholH_mat_v.matrix, &H_mat_v.matrix);  	
  /* Calculate mean vector */
  nat_grad_step(&x_vec_v.vector, &cholH_mat_v.matrix, &dfx_vec_v.vector, &mean_vec_v.vector,stepsize);

  /* Random sample from multivariate normal with 
     mean mean_vec and precision H_mat_v */
  gsl_matrix_scale(&cholH_mat_v.matrix, 1.0/sqrt(stepsize) );
  mnv_norm_rnd_cholPr(rng, &mean_vec_v.vector, &cholH_mat_v.matrix, &new_x_v.vector);

  /* Calculate propbability of new state given old state */
  double pNgO = log_mv_norm_pdf_cholP(&new_x_v.vector, &mean_vec_v.vector, &cholH_mat_v.matrix);
  
  /* evaluate model and get new state */
  double new_fx;
  int res = model->Likelihood(state->new_x,
			      model->m_params,
			      &new_fx,
			      state->new_dfx,
			      state->new_Hfx);
  if (res != GSL_SUCCESS){
    fprintf(stderr,"[smmala warning]: Model cannot be evaluated with this argument: LogLikelihood=-inf\n");
    *acc = 0;
  } else {
    gsl_vector_view new_dfx_v = gsl_vector_view_array(state->new_dfx, n);
    gsl_matrix_view new_Hfx_v = gsl_matrix_view_array(state->new_Hfx, n, n);
    
    
    res = gsl_linalg_cholesky_decomp(&new_Hfx_v.matrix);
    if (res != 0){
	fprintf(stderr,"Warning: matrix not positive definite in smmala_sample.\n Will replace matrix with I (gsl_matrix_set_identity).\n");
	gsl_matrix_set_identity (&new_Hfx_v.matrix);
    }
    gsl_matrix_memcpy(&cholH_mat_v.matrix, &new_Hfx_v.matrix);	
    /* Calculate new mean */
    nat_grad_step(&new_x_v.vector, &cholH_mat_v.matrix, &new_dfx_v.vector, &mean_vec_v.vector, stepsize);
    
    /* Calculate propbability of old state given new state */
    gsl_matrix_scale(&cholH_mat_v.matrix, 1.0/sqrt(stepsize) );
    double pOgN = log_mv_norm_pdf_cholP(&x_vec_v.vector, &mean_vec_v.vector, &cholH_mat_v.matrix);
    
    /* Accept/Reject new state */
    double mh_ratio = new_fx + pOgN - state->fx - pNgO;
    double rand_dec = log(gsl_rng_uniform(rng));
    if ( (mh_ratio > 0.0)||(mh_ratio > rand_dec) ) {
      *acc = 1;
      double* tmp;
        SWAP(kernel->x, state->new_x, tmp)
	SWAP(state->dfx, state->new_dfx, tmp)
	SWAP(state->Hfx, state->new_Hfx, tmp)
      state->fx  = new_fx;
    }else{
      *acc = 0;
    }
  }
  return GSL_SUCCESS;
}

static void smmala_kernel_adapt(mcmc_kernel* kernel, double acc_rate){
  smmala_params* params = (smmala_params*) kernel->kernel_params;
  double a=params->target_acceptance;
  if (acc_rate > 1.1*a) {
    params->stepsize = 1.3*params->stepsize;
  }
  else if(acc_rate < 0.9*a){
    params->stepsize = 0.8*params->stepsize;
  }  
}

static void smmala_kernel_free(mcmc_kernel* kernel){
	smmala_params* params = (smmala_params*)kernel->kernel_params;
	smmala_params_free(params);
	free(kernel->x);
	gsl_rng* rng = (gsl_rng*)kernel->rng;
	gsl_rng_free(rng);
	free(kernel);
}

static void smmala_kernel_print_stats(mcmc_kernel* kernel, FILE* outStream){
	smmala_params* params = (smmala_params*)kernel->kernel_params;
	fprintf(outStream, "Step-size: %5.4g\n",params->stepsize);
}

smmala_model* smmala_model_alloc(fptrPosterior_smmala Lx, fptrPrior_rnd Prx, 
								 void* model_params){
	
	smmala_model* model = (smmala_model*) malloc( sizeof (smmala_model) );
	if (model == 0){
		/* TODO: write a proper error handler */
		fprintf(stderr,"malloc failed to allocate memory for smmala_model\n");
		return 0;
	}

	model->Likelihood	= Lx;
	model->Prior_rnd = Prx;
	model->m_params = model_params;

	return model;
}

void smmala_model_free(smmala_model* model){
	free(model);
}

mcmc_kernel* smmala_kernel_alloc(int N, double step_size, smmala_model* model_function, unsigned long int seed, double target_acceptance){

	
  smmala_params* params = smmala_params_alloc(N, step_size,target_acceptance);
  if( params == NULL ){
    /* TODO: write a propper error handler */
		fprintf(stderr,"malloc failed to allocate memory for params in smmala_kernel_alloc \n");
		return 0;
  }
  
  mcmc_kernel* kernel = (mcmc_kernel*) malloc( sizeof(mcmc_kernel) );
  if (kernel == NULL){
    /* TODO: write a proper error handler */
    smmala_params_free(params);
    fprintf(stderr,"malloc failed to allocate memory for mcmc_kernel in smmala_alloc.\n");
    return 0;
  }
  
  kernel->x = (double*) malloc( N * sizeof(double) );
  if (kernel->x == NULL){
    /* TODO: write a proper error handler */
    smmala_params_free(params);
    free(kernel);
    fprintf(stderr,"malloc failed to allocate memory for mcnc_kernel.x in smmala_alloc \n");
    return 0;
  }
  
  gsl_rng_env_setup();
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  if (r == NULL){
    /* TODO: write a proper error handler */
    free(kernel->x);
    smmala_params_free(params);
    free(kernel);
    fprintf(stderr,"failed to create gsl_rng in smmala_alloc \n");
    return 0;
  }
  
  gsl_rng_set(r, seed);
  kernel->rng = r;
  
  kernel->N = N;
  kernel->model_function = model_function;
  kernel->kernel_params = params;
  kernel->ExchangeInformation = &smmala_exchange_information;
  kernel->SwapChains = &smmala_swap_chains;
  kernel->Sample = &smmala_kernel_sample;
  kernel->Adapt = &smmala_kernel_adapt;
  kernel->Init = &smmala_kernel_init;
  kernel->InitR = &smmala_kernel_init_rand;
  kernel->Free = &smmala_kernel_free;
  kernel->PrintStats = &smmala_kernel_print_stats;
  return kernel;
}
