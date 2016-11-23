#include <string.h>
#include "read_cnf.h"
//#include "model_parameters_smmala.h"

int count_rc(FILE *cnf, const regex_t *end, const regex_t *comment, int *rows, int *columns){
  int bsize=1024;
  char buffer[bsize];
  char *c;
  size_t nm=3;
  regmatch_t match[nm];
  char *fgets_return;
  rows=0;
  columns=0;
  // skip empty lines as well? 
  // skip all comment lines:
  do fgets_return=fgets(buffer,bsize,cnf);
  while (regexec(&comment[0],buffer,0,NULL,0)==0 && fgets_return!=NULL);
  // terminate line at comment symbol, if content and comments are mixed;
  //          this: "12 13 10 # these are very precise\n"
  // transforms to: "12 13 10 \0" for column counting.
  if (regexec(&comment[1],buffer,nm,match,0)==0) buffer[match[1].rm_so]='\0';
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
      buffer[match[1].rm_so]='\0';
    }
    // non empty lines which don't appear to be comments are counted as rows:
    if (regexec(&comment[0],buffer,0,NULL,0)==REG_NOMATCH) rows[0]++;
  } while (regexec(end,buffer,0,NULL,0)==REG_NOMATCH);
  //  printf("# counted %i rows in block.\n",l);
  return GSL_SUCCESS;
}

int read_block(int rows, int columns, FILE *f, int target_type, void *target, const regex_t *comment){
  int r,c;
  int bsize=1024;
  char buffer[bsize];
  size_t nm=3;
  regmatch_t match[nm];
  int match_err;
  int *v_i=NULL;
  double *v_d=NULL;
  char *line;

  if (target_type==DOUBLE_BLOCK){
    v_d=(double *) target;
  } else if (target_type==INTEGER_BLOCK){
    v_i=(int *) target;
  } else {
    fprintf(stderr,"unexpected target type in «read_block»: %i\n",target_type);
    exit(-3);
  }
  
  for (r=0;r<rows;r++){ 
    do{
      fgets(buffer,bsize,f);
      //printf("%s (is comment: %i)\n",buffer,regexec(comment,buffer,0,NULL,0)==0);
    } while (regexec(&comment[0],buffer,0,NULL,0)==0 || strlen(buffer)<2);
    //printf("⟨");
    match_err=regexec(&comment[1],buffer,nm,match,0);
    // terminate string at beginning of comment:
    if (match_err==0) buffer[match[1].rm_so]='\0';
    line=buffer;
    if (target_type==DOUBLE_BLOCK && v_d!=NULL){
      for (c=0; c<columns; c++) v_d[r*columns+c]=strtod(line,&line);
    } else if (target_type==INTEGER_BLOCK && v_i!=NULL){
      for (c=0; c<columns; c++) v_i[r*columns+c]=strtol(line,&line,10);    
    } else {
      fprintf(stderr,"read_block error for target type: %i\n",target_type);
      exit(-3);
    }
    //printf("⟩");
  }
  return EXIT_SUCCESS;
}


/* 
//just in case
int gsl_vector_square(gsl_vector *v){
  //power must be 2 or 
  int i,n=v->size;
  double s;
  for (i=0;i<n;i++) {
    s=gsl_vector_get(v,i);
    gsl_vector_set(v,i,s*s);
  }
  return GSL_SUCCESS;
}

int gsl_vector_sqrt(gsl_vector *v){
  //power must be 2 or 
  int i,n=v->size;
  double s;
  for (i=0;i<n;i++) {
    s=gsl_vector_get(v,i);
    gsl_vector_set(v,i,sqrt(s));
  }
  return GSL_SUCCESS;
}
*/

/* All normalisation functions in this file use absolte values of
 * derivatives and add them weighted by individual sds. This works
 * under the assumption that the result's standard deviation is more
 * accurately estimated by: sd(r=d/s) = |r/d|*sd(d) + |r/s|*sd(s).
 * This overestimates uncorrelated errors slightly compared to
 * (sqrt((r/d)²*sd(d))² + (r/s)²*sd(s)²)), which is the Taylor series based
 * expression.  However, in systems biology, the measurement errors
 * are often large, so the bias in a truncated Taylor series is
 * significant. Therefore, we err on the side of caution.  We also assume
 * both d and s to be positive, hence absolute values by default.
 */

int normalise_by_timepoint_with_sd(ode_model_parameters *omp, problem_size *ps){
  /* here normalisation is just one entry, containing the
   * normalisation time index. Many of the operation are done in
   * place, so variable interpretation changes.
   */
  int c,j,l;
  int C=ps->C;
  int T=ps->T;
  int F=ps->F;
  gsl_vector *d,*sd_d,*sd_s,*s;
  gsl_vector *tmp; // intermediate results
  tmp=gsl_vector_alloc(F);
  
  
  // data vector at time t_j and under condition u_c: data[c*T+j] 
  for (c=0;c<C;c++){
    l=gsl_matrix_int_get(omp->norm_t,c,0);
    for (j=0;j<T;j++) {
      d=omp->data[c*T+j];
      sd_d=omp->sd_data[c*T+j];
      s=omp->data[c*T+l];
      sd_s=omp->sd_data[c*T+l];
      // data will be scaled by 1/s
      gsl_vector_div(d,s);
      // calculate standard deviation of result:
      gsl_vector_set_zero(tmp);
      gsl_vector_add(tmp,d);
      // add both error terms
      gsl_vector_div(sd_s,s);
      gsl_vector_mul(tmp,sd_s); // so: (d/s) * (sd_s/s) = (d/s²) sd_s   
        
      gsl_vector_div(sd_d,s);
      gsl_vector_add(sd_d,tmp); // so: (sd_d/s) + sd_s × (d/s²) = (d/s)×(sd_d/d) + (d/s)×(sd_s/s);
    }
  }
  /* all datapoints are now normalised in place. Some data points are
   * now exactly 1, but the simulations will be as well, so they won't
   * contribute to the likelihood.
   */
  gsl_vector_free(tmp);
  return GSL_SUCCESS;  
}

int normalise_by_state_var_with_sd(ode_model_parameters *omp, problem_size *ps){
  /* here normalisation consists of two lines. Line 1 lists the state
   * variables to use for normalisation. Line 2 selects the time index
   * of that state variable to normalise at.
   */
  int c,i,j;
  gsl_vector *tmp; // intermediate results
  //gsl_vector_view D,SD;
  gsl_vector *s,*sd_s,*d,*sd_d;
  int i_f, i_t;
  int C=ps->C;
  int T=ps->T;
  int F=ps->F;
  // here, we need to allocate some memory for s and sd_s, because the
  // normalisation is more complex and different from state variable
  // to state variable.
  tmp=gsl_vector_alloc(F);
  s=gsl_vector_alloc(F);
  sd_s=gsl_vector_alloc(F);
  for (c=0;c<C;c++) {
    // s is defined per block of experimental conditions (c=const.)
    for (i=0;i<F;i++) {
      // fill s with the appropriate data point:
      // that is, state variable with index i_d
      i_f=gsl_matrix_int_get(omp->norm_f,c,i); 
      // at time index i_t, as specified by the normalisation matrix
      i_t=gsl_matrix_int_get(omp->norm_t,c,i);
      d=omp->data[c*T+i_t];
      sd_d=omp->data[c*T+i_t];
      gsl_vector_set(s,i,gsl_vector_get(d,i_f));
      gsl_vector_set(sd_s,i,gsl_vector_get(sd_d,i_f));
    }
    for (j=0;j<T;j++) { 
      d=omp->data[c*T+j];
      sd_d=omp->data[c*T+j];
      // data will be scaled by 1/s
      gsl_vector_div(d,s);
      // calculate standard deviation of result:
      gsl_vector_set_zero(tmp);
      gsl_vector_add(tmp,d);
      // add both error terms
      gsl_vector_div(sd_s,s);
      gsl_vector_mul(tmp,sd_s);
      // so: (d/s) * (sd_s/s) = (d/s²) sd_s   
      gsl_vector_div(sd_d,s);
      gsl_vector_add(sd_d,tmp);
      // so: (sd_d/s) + sd_s × (d/s²) = (d/s)×(sd_d/d) + (d/s)×(sd_s/s);
    }
  }
  gsl_vector_free(tmp);
  gsl_vector_free(s);
  gsl_vector_free(sd_s);
  return GSL_SUCCESS;  
}


int ratio_with_sd(gsl_matrix_sd *A, gsl_matrix_sd *B){
  /* we take the ratio A./B (by elements) B is repeated to match the
   * size of A. Both A and B contain a struct member sd, its standard
   * deviation values. A must have the same number of columns as B and
   * an integer multiple of rows, e.g. A=[B;B;B] (GNU Octave
   * Notation).  the result is stored in A.  All matrices are expected
   * to contain only positive elements.
   */
  int l;
  int C,T,F;
  gsl_matrix_view a;
  gsl_matrix *R; // for the relative error of B

  T=B->M->size1;
  C=(A->M->size1)/T;
  F=A->M->size2;
  //  printf("sizes: T=%i, C=%i, F=%i\n",T,C,F);
  R=gsl_matrix_alloc(C*T,F);
  
  
  for (l=0;l<C;l++){
    a=gsl_matrix_submatrix(A->M,l*T,0,T,F);
    gsl_matrix_div_elements(&(a.matrix),B->M); //A'=A./B
  }
  
  for (l=0;l<C;l++){
    a=gsl_matrix_submatrix(A->sd,l*T,0,T,F);//sd(A)'=sd(A)./B
    gsl_matrix_div_elements(&(a.matrix),B->M);
  }
  
  gsl_matrix_div_elements(B->sd,B->M); //sd(B)'=sd(B)./B
  gsl_matrix_memcpy(R,A->M); // copy: R=A'=A./B
  for (l=0;l<C;l++){
    a=gsl_matrix_submatrix(R,l*T,0,T,F);
    gsl_matrix_mul_elements(&(a.matrix),B->sd); //R effectively ends up being A./B*sdB/B
  }
  gsl_matrix_add(A->sd,R);
  gsl_matrix_free(R);
  return EXIT_SUCCESS;
}

int determine_problem_size(FILE *cnf, const field_expression *fe, const regex_t *comment, problem_size *ps, main_options *cnf_options){
  /* Counts columns and rows of supplied data blocks to allocate the
   * right amount of memory later on. The function returns an
   * indicator for the normalisation type.
   */
  int D=0,U=0,C=0,F=0,T=0,N=0,R=0;
  int n_f_[2]={0,0};
  int n_t_[2]={0,0};
  int rm,l;
  
  int bsize=2048;
  char buffer[bsize];
  char *var_value, *var_value_newline;
  const field_expression *current=fe;
  size_t nm=3;
  regmatch_t match[nm];
  int normalisation_type=DATA_IS_ABSOLUTE;
  
  while (!feof(cnf)){
    do {
      fgets(buffer,bsize,cnf);
    } while (regexec(&comment[0],buffer,0,NULL,0)==0);
    if (regexec(&comment[1],buffer,nm,match,0)==0) buffer[match[1].rm_so]='\0';    
    // if the string looks like this: [var_name]=[var_value]\n
    // then this will return a pointer to "="
    var_value=strchr(buffer,'='); 
    if (var_value!=NULL) { // we have a variable definition
      var_value[0]='\0'; //mark the end of the variable name
      /*printf("var_name=%s\n",buffer);
        printf("var_value=%s\n",var_value+1);
       *printf("removing newline....\n");
       */
      var_value_newline=strchr(++var_value,'\n');
      if (var_value_newline!=NULL) var_value_newline[0]='\0';      
      printf("# [cfg] {%s} ← {%s}\n",buffer,var_value);
      // set some options that can also be set from command line.
      // there is no flag whether an option was already set on command
      // line (to do later?)  at the moment. Values are compared
      // against their defaults, if an option is still at default it
      // will be changed here. If it isn't, then it must have
      // been changed on command line and the value in the file will
      // be disregarded. Command line options have higher precedence.
      if (strcmp(cnf_options->output_file,"sample.dat")==0 && strcmp(buffer,"output")==0) strcpy(cnf_options->output_file,var_value);
      else if (cnf_options->sample_size<0 && strcmp(buffer,"sample_size")==0) cnf_options->sample_size=strtol(var_value,NULL,0);
      else if (cnf_options->target_acceptance<0 && strcmp(buffer,"acceptance")==0) cnf_options->target_acceptance=strtod(var_value,NULL);
      else if (cnf_options->initial_stepsize<0 && strcmp(buffer,"step_size")==0) cnf_options->initial_stepsize=strtod(var_value,NULL);
      else if (strcmp(buffer,"t0")==0) {cnf_options->t0=strtod(var_value,NULL); printf("# t0 = %f\n",cnf_options->t0);}
    }
    else { 
      while (current != NULL){
	if (regexec(current->opening_bracket,buffer,0,NULL,0)==0){
	  /* printf("found match with regular expression %i.\n",i); */
	  /* printf(buffer); */
	  switch (current->id){ // counting columns and/or rows
	  case i_time:
	    count_rc(cnf, current->closing_bracket, comment, &T, &rm);
	    break;
	  case i_reference_input:
	    normalisation_type=DATA_NORMALISED_BY_REFERENCE;
	    count_rc(cnf, current->closing_bracket, comment, &R, &U);
	    printf("data is normalised by reference data set.\n");
	    break;
	  case i_input:
	    count_rc(cnf, current->closing_bracket, comment, &C, &U);
	    // printf("# input field has %i lines. (%i experimental conditions)\n",l,C);
	    break;
	  case i_data:
	    count_rc(cnf, current->closing_bracket, comment, &l, &F);
	    printf("# data has %i (%i × %i) lines.\n",l,C,T);
	    break;
	  case i_prior_mu:
	    count_rc(cnf, current->closing_bracket, comment, &D, &rm);
	    break;
	  case i_output:
	    count_rc(cnf, current->closing_bracket, comment, &F, &N);
	    break;
	  case i_norm_f:
	    count_rc(cnf, current->closing_bracket, comment, &n_f_[0], &n_f_[1]);
	    break;
	  case i_norm_t:
	    count_rc(cnf, current->closing_bracket, comment, &n_t_[0], &n_t_[1]);
	    break;
	    
	  }// switch
	}// if match
	current=current->next;
      }// while (regular expressions from stack)
    } //if [key]=[value]
  }// while !EOF
  printf("# file read once to determine the size of the problem.\n");
  printf("# D=%i\tN=%i\tF=%i\tU=%i\tC=%i\tR=%i\tT=%i\n",D,N,F,U,C,R,T);

  if (normalisation_type!=DATA_NORMALISED_BY_REFERENCE){
    if (n_f_[0] == 0 && n_t_[0] == C) {
      normalisation_type=DATA_NORMALISED_BY_TIMEPOINT;
      if (n_t_[0] == 1)
	R=1;
      else if (n_t_[0] == C && n_t_[1] ==1){
	R=C;
      } else {
	fprintf(stderr,"wrong number of normalisation time indices: (%i×%i)\nMust be either (1×1) or (C×1).\n",n_t_[0],n_t_[1]);
      }
    } else {
      normalisation_type=DATA_NORMALISED_BY_STATE_VAR;
      if (n_f_[1]==F && n_t_[1]==F){
	if (n_f_[0]==1 && n_t_[0]==1) {
	  R=1;
	} else if (n_f_[0]==C && n_t_[0]==C) {
	  R=C;
	} else {
	  fprintf(stderr,"Both norm_f and norm_t must have 1 or C rows: norm_f has %i rows, norm_t has %i rows.\n",n_f_[0],n_t_[0]);
	}
      }else{
	fprintf(stderr,"norm_f and norm_t must have the same number of columns, F: norm_f has %i columns, norm_t has %i columns.\n",n_f_[1],n_t_[1]);
      }
    }
  }

  //save the problem size determined from configuration file:
  ps->D=D;
  ps->C=C;
  ps->U=U;
  ps->T=T;
  // these are known from the model file:
  N=ps->N;
  // this can be checked for consistency:
  if (R==1 || R==C){
    ps->R=R;
  } else {
    fprintf(stderr,"the number of different (%i×%i) reference measurement sets must be either 1 or %i.\n",T,F,C);
    exit(-3);
  }
  if (F!=ps->F) {
    fprintf(stderr,"data file has a different number of outputs than model.\n");
    exit(-1);
  }
  if (ps->P!=D+U){
    fprintf(stderr,"error: P!=D+U.\n The number of model parameters P is not equal to\n the number of unknown paramerts D plus number of input parameters U.\n");
    exit(-2);
  }
  return normalisation_type;
}


int read_problem_definition(FILE *cnf, ode_model_parameters *omp, gsl_matrix_sd *RD, const field_expression *fe, const regex_t *comment, problem_size *ps, main_options *cnf_options){
  int bsize=2048;
  char buffer[bsize];
  char *c;
  //  size_t nm=3;
  //  regmatch_t match[nm];
  const field_expression *current;
  while (!feof(cnf)){
    // discard all comment lines
    do c=fgets(buffer,bsize,cnf);
    while (regexec(&comment[0],buffer,0,NULL,0)==0 && c!=NULL);
    /* now that the line is not a comment or empty:
     * reset the stack pointer to the top and try all field
     * expressions on that line
     */
    current=fe;
    while (current!=NULL){
      if (regexec(current->opening_bracket,buffer,0,NULL,0)==0){
        /* printf("found match with regular expression id %i.\n",i); */
        /* printf(buffer); */
	switch (current->id){// reading the data
	case i_time:
	  read_block(ps->T,1,cnf,DOUBLE_BLOCK,omp->t->data,comment);
	  printf("# measurement time(s) read.\n");
	  break;
	case i_reference_input:
	  //printf("# reading reference input.\n");
	  read_block(1,ps->U,cnf,DOUBLE_BLOCK,omp->reference_u->data,comment);
	  printf("# reference input read.\n");
	  break;
	case i_reference_data:
	  read_block(ps->T,ps->F,cnf,DOUBLE_BLOCK,RD->M->data,comment);
	  printf("# reference data read.\n");
	  break;
	case i_sd_reference_data:
	  read_block(ps->T,ps->F,cnf,DOUBLE_BLOCK,RD->sd->data,comment);
	  printf("# standard deviation of reference data read.\n");
	  break;
	case i_input:
	  read_block(ps->C,ps->U,cnf,DOUBLE_BLOCK,omp->input_u->data,comment);
	  printf("# input read.\n");
	  break;
	case i_data:
	  read_block((ps->C)*(ps->T),ps->F,cnf,DOUBLE_BLOCK,omp->Data->data,comment);
	  printf("# data read.\n");
	  break;
	case i_sd_data:
	  read_block((ps->C)*(ps->T),ps->F,cnf,DOUBLE_BLOCK,omp->sdData->data,comment);
	  printf("# standard deviation of data read.\n");
	  break;
	case i_prior_mu:
	  read_block(ps->D,1,cnf,DOUBLE_BLOCK,omp->prior_mu->data,comment);
	  printf("# prior mean read.\n");
	  break;
	case i_prior_icov:
	  read_block(ps->D,ps->D,cnf,DOUBLE_BLOCK,omp->prior_inverse_cov->data,comment);
	  printf("# prior inverse covariance matrix read.\n");
	  break;
	case i_output:
	  read_block(ps->F,ps->N,cnf,DOUBLE_BLOCK,omp->output_C->data,comment);
	  break;
	case i_norm_f:
	  read_block(ps->R,ps->F,cnf,INTEGER_BLOCK,omp->norm_f,comment);
	  break;
	case i_norm_t:
	  read_block(ps->R,ps->F,cnf,INTEGER_BLOCK,omp->norm_t,comment);
	  break;
	}// switch what to do.
	break; // if matched: break, else try the next regular expression.
      } else {
	current=current->next;
      }// if match
    }// while (regular expressions)
  }// while !EOF
  return EXIT_SUCCESS;
}


field_expression* field_expression_stack(int id, field_expression *top, regex_t *open, regex_t *close){
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
  fe->id=id;
  fe->next=top; //old top
  // this item is now the new top of the stack;
  return fe;
}

int field_names_init(field_names *fn){
  // field names should not be longer than 128 characters.
  // n fields are known to the parser
  int i; 
  fn->n=12;
  fn->max_length=72;
  fn->name=(char **) calloc(fn->n,sizeof(char*));
  for (i=0;i<fn->n;i++) {
    fn->name[i]=(char*) calloc(fn->max_length,sizeof(char));
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
  strcpy(fn->name[i_norm_f],"norm_f");
  strcpy(fn->name[i_norm_t],"norm_t");

  return EXIT_SUCCESS;
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
  fe=NULL;
  for (i=0;i<n;i++){
    // we assume that fe is completely uninitialised
    opening = (regex_t*) malloc(sizeof(regex_t));
    closing = (regex_t*) malloc(sizeof(regex_t));
    //regcomp(opening,"\\[time\\]",REG_EXTENDED);
    //regcomp(closing,"\\[/time\\]",REG_EXTENDED);
    sprintf(sptr,"\\[%s\\]",fn->name[i]);
    regcomp(opening,sptr,REG_EXTENDED);
    sprintf(sptr,"\\[/%s\\]",fn->name[i]);
    regcomp(closing,sptr,REG_EXTENDED);  
    fe=field_expression_stack(i,fe,opening,closing);
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
  //  int D=0,U=0,C=0,F=0,T=0,N=0;
  int c,i;
  regex_t comment[2];
  field_expression *fe;
  field_names fn;

  field_names_init(&fn);
  /*comments in file: a line that is only a comment: "_blank_ #
   * comment" must be disregarded when line counting; a line that
   * contains both content and comments must be read up to the comment
   * symbol.
   */
  regcomp(&comment[0],"^\\s*((#|//|%)|$)",REG_EXTENDED);
  regcomp(&comment[1],"\\w+\\s*(#|//|%)",REG_EXTENDED);
  
  fe=field_expression_init(&fn);

  omp->normalisation_type=determine_problem_size(cnf,fe, comment, ps, cnf_options);
  /* now we must allocate the necessary memory and fill it with values
   */
  ode_model_parameters_alloc(omp, ps);
  omp->t0=cnf_options->t0;
  omp->D=ps->D;
  printf("# smmala memory allocated.\n");

  // reference data:
  gsl_matrix_sd RD;
  RD.M=gsl_matrix_alloc(ps->T,ps->F);
  RD.sd=gsl_matrix_alloc(ps->T,ps->F);
  gsl_matrix_sd ND; // Data to be Normalised by reference
  ND.M=omp->Data;
  ND.sd=omp->sdData;
  int R=ps->R;
  int C=ps->C;
  rewind(cnf);
  printf("# rewind.\n");
  read_problem_definition(cnf, omp, &RD, fe, comment, ps, cnf_options);
  
  printf("# configuration read.\n");

  /* now we perform the data normalisation
   */
  printf("# calculating relative data.");

  switch (omp->normalisation_type){
  case DATA_NORMALISED_BY_REFERENCE:
    printf("# taking ratios of Data and Refrence Data...");
    ratio_with_sd(&ND,&RD); //Data/Reference Data.
    printf("# done.\n");
    break;
  case DATA_NORMALISED_BY_TIMEPOINT:
    printf("# normalising using one time instance...");
    if (R<C){ // user supplied fewer than C normalisation time indices
      for (c=R;c<C;c++) { //copy this time index for all experimental conditions 
	gsl_matrix_int_set(omp->norm_t,c,0,gsl_matrix_int_get(omp->norm_t,c%R,0));
      }
    }
    normalise_by_timepoint_with_sd(omp,ps);
    printf("# done.\n");
    break;
  case DATA_NORMALISED_BY_STATE_VAR:
    printf("# normalising using another state variable at specified time instance...");
    if (R<C){ // user supplied just one normalisation line
      for (c=R;c<ps->C;c++) { //copy this line for all c 
	for (i=0;i<ps->F;i++){
	  gsl_matrix_int_set(omp->norm_t,c,i,gsl_matrix_int_get(omp->norm_t,c%R,i));
	  gsl_matrix_int_set(omp->norm_f,c,i,gsl_matrix_int_get(omp->norm_f,c%R,i));
	}
      }
    }
    normalise_by_state_var_with_sd(omp,ps);
    printf("# done.\n");
    break;
  }
  // cleanup
  gsl_matrix_free(RD.M);
  gsl_matrix_free(RD.sd);
  while (fe!=NULL) { //empty stack
    regfree(fe->opening_bracket);
    regfree(fe->closing_bracket);
    fe=fe->next;    
  }

  return EXIT_SUCCESS;
}
