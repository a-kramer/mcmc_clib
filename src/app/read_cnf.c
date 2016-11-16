#include <string.h>
#include "read_cnf.h"
#include "model_parameters_smmala.h"

#define i_time 0
#define i_reference_input 1
#define i_reference_data 2
#define i_sd_reference_data 3
#define i_input 4
#define i_data 5
#define i_sd_data 6
#define i_prior_mu 7
#define i_prior_icov 8
#define i_output 9
#define i_norm 10

int count_rc(FILE *cnf, regex_t *end, regex_t *comment, int *rows, int *columns){
  int bsize=1024;
  char buffer[bsize];
  char *c;
  size_t nm=3;
  regmatch_t match[nm];

  rows=0;
  columns=0;
  // skip empty lines as well? 
  // skip all comment lines:
  do fgets(buffer,bsize,cnf);
  while (regexec(&comment[0],buffer,0,NULL,0)==0);
  // terminate line at comment symbol, if content and comments are mixed;
  //          this: "12 13 10 # these are very precise\n"
  // transforms to: "12 13 10 \0" for column counting.
  if (regexec(&comment[1],buffer,nm,match,0)==0) buffer[match[1].so]='\0';
  c=buffer;

  printf("counting columns of «%s»...",c);
  while ( memchr("\n\0",c[0],2)==NULL ){
    /* move forward until we find a word character (not space or tab) */
    while ( isspace(c[0]) ) c++;
    columns[0]++; // now that we found one, we increase the column counter;
    while ( isalnum(c[0]) || ispunct(c[0])) c++;
  };

  do { 
    if (fgets(buffer,bsize,cnf)==NULL) {
      printf("error reading config file: end of file reached; missing closing tag?\n");
      exit(1);
    } else if ( regexec(&comment[1],buffer,nm,match,0)==0) {
      buffer[match[1].so]='\0';
    }
    // non empty lines which don't appear to be comments are counted as rows:
    if (regexec(&comment[0],buffer,0,NULL,0)==REG_NOMATCH) rows[0]++;
  } while (regexec(end,buffer,0,NULL,0)==REG_NOMATCH);
  //  printf("# counted %i rows in block.\n",l);
  return GSL_SUCCESS;
}

int read_vector(char *c, double *vector, const int length){
  int l=0;
  //printf("read_");
  for (l=0; l<length; l++) vector[l]=strtod(c,&c);
  //printf("columns");
  return l;
}

int read_block(int rows, int columns, FILE *f, double *target, regex_t **comment){
  int l;
  int bsize=1024;
  char buffer[bsize];
  size_t nm=3;
  regmatch_t match[nm];
  int match_err;
  for (l=0;l<rows;l++){ 
    do{
      fgets(buffer,bsize,f);
      //printf("%s (is comment: %i)\n",buffer,regexec(comment,buffer,0,NULL,0)==0 || strlen(buffer)<2);
    } while (regexec(comment[0],buffer,0,NULL,0)==0 || strlen(buffer)<2);
    //printf("⟨");
    match_err=regexec(comment[1],buffer,nm,match,0);
    // terminate string at beginning of comment:
    if (match_err==0) buffer[match[1].so]='\0';
    read_vector(buffer,&target[l*columns],columns);
    //printf("⟩");
  }
  return EXIT_SUCCESS;
}

int ratio_with_sd(gsl_matrix *A, gsl_matrix *sdA, gsl_matrix *B, gsl_matrix *sdB){
  /* we take the ratio A./B (by elements) B is repeated to match the
   * size of A. A must have the same number of columns as B and an
   * integer multiple of rows, e.g. A=[B;B;B] (GNU Octave Notation).
   * the result is stored in A.  All matrices are expected to contain
   * only positive elements.
   */
  int l;
  int C,T,F;
  gsl_matrix_view a;
  gsl_matrix *R; // for the relative error of B

  T=B->size1;
  C=(A->size1)/T;
  F=A->size2;
  //  printf("sizes: T=%i, C=%i, F=%i\n",T,C,F);
  R=gsl_matrix_alloc(C*T,F);
  
  
  for (l=0;l<C;l++){
    a=gsl_matrix_submatrix(A,l*T,0,T,F);
    gsl_matrix_div_elements(&(a.matrix),B); //A'=A./B
  }
  
  for (l=0;l<C;l++){
    a=gsl_matrix_submatrix(sdA,l*T,0,T,F);//sdA'=sdA./B
    gsl_matrix_div_elements(&(a.matrix),B);
  }
  
  gsl_matrix_div_elements(sdB,B); //sdB'=sdB./B
  gsl_matrix_memcpy(R,A); // copy: R=A'=A./B
  for (l=0;l<C;l++){
    a=gsl_matrix_submatrix(R,l*T,0,T,F);
    gsl_matrix_mul_elements(&(a.matrix),sdB); //R effectively ends up being A./B*sdB/B
  }
  gsl_matrix_add(sdA,R);
  gsl_matrix_free(R);
  return EXIT_SUCCESS;
}

int determine_problem_size(FILE *cnf, const field_expressions *fe, const regex_t *comment, problem_size *ps, main_options *cnf_options){
  // file content buffer:
  int D=0,U=0,C=0,F=0,T=0,N=0;
  int i,rm,l;
  int n=fe->n; // number of fields
  regex_t *cpt=fe->cpt;
  int bsize=2048;
  char buffer[bsize];
  char *var_value, *var_value_newline;
  field_expression *current=fe;
  size_t nm=3;
  regmatch_t match[nm];

  while (!feof(cnf)){
    do {
      fgets(buffer,bsize,cnf);
    } while (regexec(&comment[0],buffer,0,NULL,0)==0);
    if (regexec(&comment[1],buffer,nm,match,0)==0) buffer[match[1].so]='\0';    
    var_value=strchr(buffer,'='); // string looks like this: [var_name]=[var_value]
    if (var_value!=NULL) { // we have a variable definition
      var_value[0]='\0'; //mark the end of the variable name
      /*printf("var_name=%s\n",buffer);
        printf("var_value=%s\n",var_value+1);
       *printf("removing newline....\n");
       */
      var_value_newline=strchr(++var_value,'\n');
      if (var_value_newline!=NULL) var_value_newline[0]='\0';      
      printf("# [cfg] {%s} ← {%s}\n",buffer,var_value);
      
      if (strcmp(cnf_options->output_file,"sample.dat")==0 && strcmp(buffer,"output")==0) strcpy(cnf_options->output_file,var_value);
      else if (cnf_options->sample_size<0 && strcmp(buffer,"sample_size")==0) cnf_options->sample_size=strtol(var_value,NULL,0);
      else if (cnf_options->target_acceptance<0 && strcmp(buffer,"acceptance")==0) cnf_options->target_acceptance=strtod(var_value,NULL);
      else if (cnf_options->initial_stepsize<0 && strcmp(buffer,"step_size")==0) cnf_options->initial_stepsize=strtod(var_value,NULL);
      else if (strcmp(buffer,"t0")==0) {omp->t0=strtod(var_value,NULL); printf("# t0 = %f\n",omp->t0);}
    }
    else { 
      for (i=0;i<n;i++){
	if (regexec(fe->[i],buffer,0,NULL,0)==0){
	  /* printf("found match with regular expression %i.\n",i); */
	  /* printf(buffer); */
	  switch (i){ // counting columns and/or rows
	  case 0: /*time*/
	    count_rc(cnf, &cpt[n+i], comment, &T, &rm);
	    break;
	  case 1: /* reference input */
	    omp->data_is_relative=1;
	    printf("data is normalised by reference data set.\n");
	    break;
	  case 4: /* inputs */
	    count_rc(cnf, &cpt[n+i], comment, &C, &U);
	    //	  printf("# input field has %i lines. (%i experimental conditions)\n",l,C);
	    break;
	  case 5: /* data */
	    count_rows(cnf, &cpt[n+i], comment, &l, &F);
	    printf("# data has %i (%i × %i) lines.\n",l,C,T);
	    break;
	  case 7: /* problem size */
	    count_rc(cnf, &cpt[n+i], comment,&D,&rm);
	    break;
	  case 9: /* output structure */
	    count_rows(cnf, &cpt[n+i], comment,&F,&N);
	    break;
	  case 10: /* normalisation matrix*/
	    do fgets(buffer,bsize,cnf);
	    while (regexec(comment[0],buffer,0,NULL,0)==0);
	    if (regexec(comment[1],buffer,nm,match,0)==0) buffer[match[1].so]='\0';
	    ps->n2=count_columns(buffer);
	    ps->n1=count_rows(cnf,&cpt[n+i],&comment);
	  }// switch
	  break;
	}// if match
      }// for (regular expressions)
    }// while !EOF
  }
  printf("# file read once to determine the size of data and inputs.\n");
  printf("# D=%i\tN=%i\tF=%i\tU=%i\tC=%i\tT=%i\n",D,N,F,U,C,T);
  //save the problem size determined from configuration file:
  ps->D=D;
  ps->C=C;
  ps->U=U;
  ps->T=T;
  // these are known from the model file:
  P=ps->P;
  N=ps->N;
  // this can be checked for consistency:
  if (F!=ps->F) {
    fprintf(stderr,"data file has a different number of outputs than model.\n");
    return -1;
  }
  if (ps->P!=D+U){
    fprintf(stderr,"number of model parameters P is not equal to number of unknown paramerts D and number of input parameters U: P!=D+U.\n");
    return -2;
  }

}


int read_problem_definition(FILE *cnf, const field_expressions fe, const regex_t *comment, problem_size *ps, main_options *cnf_options){
  int bsize=2048;
  char buffer[bsize];
  size_t nm=3;
  regmatch_t match[nm];
  while (!feof(cnf)){
    do fgets(buffer,bsize,cnf);
    while (regexec(&comment[0],buffer,0,NULL,0)==0);
    for (i=0;i<n;i++){
      if (    regexec(&cpt[i],buffer,0,NULL,0)==0){
        /* printf("found match with regular expression %i.\n",i); */
        /* printf(buffer); */
	switch (i){// reading the data
	case 0: /*time*/
	  read_block(T,1,cnf,omp->t->data,&comment);
	  printf("# measurement time(s) read.\n");
	  break;
	case 1: /* reference input */
	  //printf("# reading reference input.\n");
	  read_block(1,U,cnf,omp->reference_u->data,&comment);
	  printf("# reference input read.\n");
	  break;
	case 2: /* reference data */
	  read_block(T,F,cnf,reference_data->data,&comment);
	  printf("# reference data read.\n");
	  break;
	case 3: /* sd reference data */
	  read_block(T,F,cnf,sd_reference_data->data,&comment);
	  printf("# standard deviation of reference data read.\n");
	  break;
	case 4: /* inputs */
	  read_block(C,U,cnf,omp->input_u->data,&comment);
	  printf("# input read.\n");
	  break;
	case 5: /* data */
	  read_block(C*T,F,cnf,omp->Data->data,&comment);
	  printf("# data read.\n");
	  break;
	case 6: /* sd data */
	  read_block(C*T,F,cnf,omp->sdData->data,&comment);
	  printf("# standard deviation of data read.\n");
	  break;
	case 7:
	  read_block(D,1,cnf,omp->prior_mu->data,&comment);
	  printf("# prior mean read.\n");
	  break;
	case 8: /* prior inverse covariance */
	  read_block(D,D,cnf,omp->prior_inverse_cov->data,&comment);
	  printf("# prior inverse covariance matrix read.\n");
	  break;
	case 9: /* output structure */
	  read_block(F,N,cnf,omp->output_C->data,&comment);
	  break;
	case 10: /* normalisation structure */
	  read_block(n1,n2,cnf,omp->normalisation,&comment);
	}// switch
	break;
      }// if match
    }// for (regular expressions)
  }// while !EOF

}


field_expression* field_expression_stack(field_expression *top, regex_t *open, regex_t *close){
  /* this creates a new struct field_expression struct with contents
   * "open" and "close". top is a pointer to the current top of the
   * stack. the new element will link to the former top of the stack.
   */
  field_expression* fe;
  // append new entries:
  fe = (field_expression*) malloc(sizeof(field_expression));
  // go to the new connected memory block;
  fe->opening_bracket=open;
  fe->closing_bracket=close;
  fe->next=top; // this is now the new end of the array;
  return fe;
}

int field_names_init(field_names *fn){
  // field names should not be longer than 128 characters.
  // n fields are known to the parser
  int i; 
  size_t n=11;
  fn->name=(char **) calloc(n,sizeof(char*));
  for (i=0;i<n;i++) {
    fn->name[i]=(char*) calloc(128,sizeof(char));
  }
  
  strcpy(fn->name[i_time],"time");
  strcpy(fn->name[i_data],"data");
  strcpy(fn->name[i_reference_data],"ref(erence)?_data");
  strcpy(fn->name[i_sd_data],"sd_data");
  strcpy(fn->name[i_sd_reference_data],"sd_ref(erence)?_data");
  strcpy(fn->name[i_input],"input");
  strcpy(fn->name[i_reference_input],"reference_input");
  strcpy(fn->name[i_prior_mu],"prior_mu");
  strcpy(fn->name[i_prior_icov],"prior_inv(erse)?_cov(ariance)?");
  strcpy(fn->name[i_output],"output");
  strcpy(fn->name[i_norm],"normali[sz]ation");
  fn->n=n;
  fn->max_length=128;
}

field_expression* field_expression_init(field_names *fn){
  //field expressions form a dynamic array with pointers between the
  //elements. The array is handled like a stack.
  // [B]->[A]->NULL.
  //new element: [C]->[B]->[A]->NULL;

  regex_t *opening, *closing;
  int i=0;
  size_t n=fn->n;
  field_expression *fe;
  char sptr[fn->max_length];

  for (i=0;i<n;i++){
    // we assume that fe is completely uninitialised
    opening = (regext_t*) malloc(sizeof(regex_t));
    closing = (regext_t*) malloc(sizeof(regex_t));
    //regcomp(opening,"\\[time\\]",REG_EXTENDED);
    //regcomp(closing,"\\[/time\\]",REG_EXTENDED);
    sprintf(sptr,"\\[%s\\]",fn->name[i]);
    regcomp(opening,sptr,REG_EXTENDED);
    sprintf(sptr,"\\[/%s\\]",fn->name[i]);
    regcomp(closing,sptr,REG_EXTENDED);  
    fe=field_expression_stack(fe,opening,closing);
  }
  
  return fe;
}

int parse_config(FILE *cnf, ode_model_parameters *omp,  problem_size *ps, main_options *cnf_options){
  /* This function parses the configuration file The configuration
   * file contains xml like expressions containing data, measurement
   * time points, etc.
   * We start by defining regular expressions for the field names
   * 
   */
  int D=0,U=0,C=0,F=0,T=0,N=0;
  int i,rm,l;
  
  int n=11; // number of fields
  // regular expression matching:
  regex_t cpt[2*n]; // defines the patterns for the field names
  regex_t comment[2];

  field_expression *fe;
  field_names fn;
  // assume data to be absolute a priori
  omp->data_is_relative=0; 

  field_names_init(&fn);
  /*comments in file: a line that is only a comment: "_blank_ #
   * comment" must be disregarded when line counting; a line that
   * contains both content and comments must be read up to the comment
   * symbol.
   */
  regcomp(&comment[0],"^\s*((#|//|%)|$)",REG_EXTENDED);
  regcomp(&comment[1],"\w+\s*(#|//|%)",REG_EXTENDED);
  
  fe=field_expression_init(&fn);

  determine_problem_size(cnf,fe, comment, ps, cnf_options);
  /* now we must allocate the necessary memory and fill it with values
   */
  gsl_matrix *reference_data;
  gsl_matrix *sd_reference_data;

  ode_model_parameters_alloc(omp, ps);

  reference_data=gsl_matrix_alloc(T,F);
  sd_reference_data=gsl_matrix_alloc(T,F);
  //  printf("# memory allocated.\n");
  /* we reset the FILE structure and parse again, now reading the
   * values.
   */
  rewind(cnf);
  printf("# rewind.\n");
  read_problem_definition(cnf,fe, comment, ps, cnf_options);
  printf("# configuration read.\n");
  /* now we take the ratios of data and reference data
   */

  if (omp->data_is_relative){
    printf("# calculating relative data (ratios)...");
    for (i=0;i<C;i++)ratio_with_sd(omp->Data[i],omp->sdData[i],reference_data,sd_reference_data);
    printf("# done.\n");
  }

  //  omp->t0=0.0;

  // cleanup
  gsl_matrix_free(reference_data);
  gsl_matrix_free(sd_reference_data);
  for (i=0;i<2*n;i++)  regfree(&cpt[i]);

  return EXIT_SUCCESS;
}
