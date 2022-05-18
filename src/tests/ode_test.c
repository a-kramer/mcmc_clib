#include <stdlib.h>
#include <stdio.h>
#include "../ode/ode_model.h"

gsl_vector* linspace(double a, double b, size_t n){
	gsl_vector *v=gsl_vector_alloc(n);
	int i;
	double d=(b-a)/(n-1);
	assert(d>0);
	assert(v && v->data);
	for (i=0;i<n;i++){
		v->data[i]=a;
		a+=d;
	}
	return v;
}

gsl_vector** gsl_vectors_alloc(size_t d, int n){
	int i;
	gsl_vector **v=malloc(sizeof(gsl_vector*)*n);
	for (i=0;i<n;i++){
		v[i]=gsl_vector_alloc(d);
	}
	return v;
}

gsl_matrix** gsl_matrices_alloc(size_t d1, size_t d2, int n){
	int i;
	gsl_matrix **m=malloc(sizeof(gsl_matrix*)*n);
	for (i=0;i<n;i++){
		m[i]=gsl_matrix_alloc(d1,d2);
	}
	return m;
}


int main(int argc, char *argv[])
{
	ode_model *M=ode_model_load_from_file("test.so");
	double t0=0.0;
	ode_ivp *ivp=ode_ivp_alloc(M);
	ode_ivp_init(ivp,t0,NULL,NULL);	

	int ny=ivp->y->size;
	int np=ivp->params->size;
	int nfy=ivp->fy->size;
	int nt=100;
	gsl_vector *t=linspace(t0,10,nt);
	gsl_vector **y=gsl_vectors_alloc(ny,nt);
	gsl_vector **fy=gsl_vectors_alloc(nfy,nt);
	gsl_matrix **yS=gsl_vectors_alloc(ny,np,nt);
	gsl_matrix **fyS=gsl_vectors_alloc(nfy,np,nt);

	int status;
	assert((status=ode_ivp_solve(ivp, t, y, fy, yS, fyS, NULL))=GSL_SUCCESS);
	/* print on screen or write to a file */
	return EXIT_SUCCESS;
}
