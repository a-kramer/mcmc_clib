#include "model_parameters_smmala.h"
#include <assert.h>
/* Data will always be stored in blocks (matrices), there will be
 * vector views to identify specific rows.	Whether omp->Data or omp->data_block
 * is used in memory allocation is up to the data read functions.
 */
int bayesian_parameters_init(bayesian_parameters *bpar, int C){
	int i;
	assert(C>0);
	bpar->E=(experiment**) malloc(sizeof(experiment*)*C);
	for (i=0;i<C;i++) {
		bpar->E[i]=(experiment*) malloc(sizeof(experiment));
		bpar->E[i]->t=NULL;
		bpar->E[i]->init_y=NULL;
		bpar->E[i]->input_u=NULL;
		bpar->E[i]->single=NULL;
	}
	bpar->Data=NULL;
	bpar->sdData=NULL;
	for (i=0;i<C;i++) {
		bpar->E[i]->data_block=NULL;
		bpar->E[i]->sd_data_block=NULL;
		bpar->E[i]->normalise=malloc(sizeof(normalisation_t));
		bpar->E[i]->view=malloc(sizeof(view_t));
	}
	bpar->S_approx=malloc(sizeof(sensitivity_approximation*)*C);
	for (i=0;i<C;i++){
		bpar->S_approx[i]=malloc(sizeof(sensitivity_approximation));
	}

	bpar->prior=malloc(sizeof(prior_t));
	//printf("# %i+1 experiment structures allocated.\n",C);
	return GSL_SUCCESS;
}


int experiment_alloc(experiment *E, size_t T){
	int j;
	//printf("[%s] P=%i.\n",__func__,P); fflush(stdout);
	E->y=(gsl_vector **) malloc(sizeof(gsl_vector*)*T);
	E->fy=(gsl_vector **) malloc(sizeof(gsl_vector*)*T);
	E->yS=(gsl_matrix **) malloc(sizeof(gsl_matrix*)*T);
	E->fyS=(gsl_matrix **) malloc(sizeof(gsl_matrix*)*T);
	E->oS=(gsl_matrix **) malloc(sizeof(gsl_matrix*)*T);

	E->normalise->fy=malloc(sizeof(gsl_vector*)*T);
	E->normalise->fyS=malloc(sizeof(gsl_matrix*)*T);
	E->normalise->data=malloc(sizeof(gsl_vector*)*T);
	E->normalise->stdv=malloc(sizeof(gsl_vector*)*T);
	return GSL_SUCCESS;
}

int bayesian_parameters_alloc(bayesian_parameters *bpar){
	/* problem_size structure contains:
	 * N: number of states
	 * D: number of model parameters
	 * F: number of output functions
	 * U: number of input parameters
	 * C: number of different experimental conditions
	 * T: number of measured time points
	 */
	int i;
	int N=bpar->size->N; // number of state variables
	int C=bpar->size->C; // number of experimental conditions
	int F=bpar->size->F; // number of output functions
	int D=bpar->size->D; // number of sampling parameters (to be estimated)
	//int U=bpar->size->U; // number of input parameters (known)
	//int T=bpar->size->T; // number of measurement time instances
	int P=bpar->size->P; // number of total parameters D+U (a consistency check
				 // between ode_model and mcmc configuration file)

	//printf("[model_parameters_alloc] allocating model parameters with: N=%i,\tP=%i,\tT=%i.\n",N,P,T);
	/* initialise gsl_vectors and matrices these will hold the
	 * differential equation data for each experimental condition c and
	 * time point t_j, like this: y[c*T+j]=gsl_vector(N)
	 */
	bpar->tmpF=gsl_vector_alloc(F);
	bpar->tmpDF=gsl_matrix_alloc(D,F);
	for(i=0;i<C;i++){
		bpar->S_approx[i]->jacobian_y=gsl_matrix_alloc(N,N);
		bpar->S_approx[i]->jacobian_p=gsl_matrix_alloc(P,N);
		bpar->S_approx[i]->R=gsl_matrix_alloc(P,N); /* working memory */
		bpar->S_approx[i]->tau=gsl_vector_alloc(N);
		bpar->S_approx[i]->x=gsl_vector_alloc(N);
		bpar->S_approx[i]->r=gsl_vector_alloc(N);
		bpar->S_approx[i]->eJt=gsl_matrix_alloc(N,N);
		bpar->S_approx[i]->Jt=gsl_matrix_alloc(N,N);
	}
	//printf("[model_parameters_alloc] y, fy,yS, fyS, oS, yS0, nfy and nfyS.\n");
	for (i=0;i<C;i++){
		experiment_alloc(bpar->E[i],bpar->size);
		bpar->E[i]->t0=bpar->t0;
	}
	bpar->p=gsl_vector_alloc(P);
	// prior
	bpar->prior->p=gsl_permutation_alloc((size_t) D);
	bpar->prior->n=3;
	bpar->prior->tmp=malloc(sizeof(gsl_vector*) * bpar->prior->n);
	for (i=0;i< bpar->prior->n;i++){
		bpar->prior->tmp[i]=gsl_vector_alloc(D);
	}
	//printf("[model_parameters_alloc] done.\n");
	return EXIT_SUCCESS;
}

int bayesian_parameters_free(bayesian_parameters *bpar){
	int i,j;
	int C=bpar->size->C;
	int T=bpar->size->T;

	for (i=0;i<C;i++){
		gsl_matrix_free(bpar->S_approx[i]->jacobian_y);
		gsl_matrix_free(bpar->S_approx[i]->jacobian_p);
		gsl_vector_free(bpar->S_approx[i]->tau);
		gsl_vector_free(bpar->S_approx[i]->x);
		gsl_vector_free(bpar->S_approx[i]->r);
		gsl_matrix_free(bpar->S_approx[i]->Jt);
		gsl_matrix_free(bpar->S_approx[i]->eJt);
		free(bpar->S_approx[i]);
	}
	free(bpar->S_approx);

	for (i=0;i<C;i++){
		free(bpar->E[i]->view->data_row);
		free(bpar->E[i]->view->sd_data_row);
		free(bpar->E[i]->data);
		free(bpar->E[i]->sd_data);
		T=bpar->E[i]->t->size;
		for (j=0;j<T;j++){
			gsl_vector_free(bpar->E[i]->y[j]);
			gsl_vector_free(bpar->E[i]->fy[j]);
			gsl_matrix_free(bpar->E[i]->yS[j]);
			gsl_matrix_free(bpar->E[i]->fyS[j]);
			gsl_matrix_free(bpar->E[i]->oS[j]);
		}
		free(bpar->E[i]);
	}
	free(bpar->E);
	gsl_vector_free(bpar->p);
	free(bpar->size);
	return EXIT_SUCCESS;
}
