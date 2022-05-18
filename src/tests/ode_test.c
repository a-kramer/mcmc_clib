#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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

void gsl_vectors_write(gsl_vector **v, size_t n, const char *f){
	assert(v && v[0]);
	FILE *F=fopen(f,"w");
	int rank = 2;
	int d=v[0]->size;
	int dim[2] = {d,n};
	int i;
	fwrite(&rank,sizeof(int),1,F);
	for (i=0;i<rank;i++){
		fwrite(&dim,sizeof(int),rank,F);
	}
	for (i=0;i<n && v[i]!=NULL;i++){
		fwrite(v[i]->data,sizeof(double),v[i]->size,F);
	}
	fclose(F);
}

void gsl_matrices_write(gsl_matrix **m, size_t n, const char *f){
	assert(m && m[0]);
	FILE *F=fopen(f,"w");
	int rank = 3;
	int dim[3] = {m[0]->size2,m[0]->size1,n};
	int i;
	fwrite(&rank,sizeof(int),1,F);
	for (i=0;i<rank;i++){
		fwrite(&dim,sizeof(int),rank,F);
	}
	for (i=0;i<n && m[i]!=NULL;i++){
		fwrite(m[i]->data,sizeof(double),m[i]->size,F);
	}
	fclose(F);
}


/* if the string s contains the character c, find returns the position
	 of c, otherwise the negated length of the string (-1*strlen(s))*/
int find(const char *s, int c)
{
	int i=0;
	assert(s);
	while (s[i]!=c && s[i]!='\0'){
		i++;
	}
	if (s[i]!=c) i*=-1;
	return i;
}

int main(int argc, char *argv[])
{
	assert(argc>1);
	int l=strlen(argv[1]);
	char *modelf=malloc(l+64);
	*(memcpy(modelf,argv[1],l+1)+l)='\0';
	ode_model *M=ode_model_load_from_file(modelf);
	assert(M);
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

	char *dot=strrchr(modelf,'.');
	if (!dot) dot=modelf+l;
	*dot='\0';

	assert((status=ode_ivp_solve(ivp, t, y, fy, yS, fyS, NULL))=GSL_SUCCESS);
	/* print on screen or write to a file */
	strcpy(dot,"_y.double");
	gsl_vectors_write(y,t->size,modelf);
	strcpy(dot,"_fy.double");
	gsl_vectors_write(fy,t->size,modelf);
	strcpy(dot,"_yS.double");
	gsl_vectors_write(yS,t->size,modelf);
	strcpy(dot,"_fyS.double");
	gsl_vectors_write(fyS,t->size,modelf);
	return EXIT_SUCCESS;
}
