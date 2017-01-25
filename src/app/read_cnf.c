#include "read_cnf.h"


gsl_object* gsl_object_alloc(){
  gsl_object *new;
  new=(gsl_object*) malloc(sizeof(gsl_object));
  new->gsl=(gsl_t*) malloc(sizeof(gsl_t));
  printf("# gsl object allocated\n");
  return new;
}

int gsl_object_init(gsl_object *G){
  G->is_double=0;
  G->is_int=0;
  G->is_matrix=0;
  G->is_vector=0;
  printf(" gsl object initialised.\n");
  return GSL_SUCCESS;
}

int gsl_object_free(gsl_object *G){
  if (G->is_matrix){
    if (G->is_double) {
      gsl_matrix_free(G->gsl->matrix);
    } else {
      gsl_matrix_int_free(G->gsl->matrix_int);
    }
  }  else {
    if (G->is_double) {
      gsl_vector_free(G->gsl->vector);
    } else {
      gsl_vector_int_free(G->gsl->vector_int);
    }
  }
  free(G->gsl);
  free(G);
  return GSL_SUCCESS;
}


var_ndim_array* read_block_nd(FILE *f, int target_type, const regex_t *end, const regex_t *comment){
  int r=0,c=0,pre_c=0; // number of rows and columns, pre_c is there
		       // to compare two subsequent rows
  int bsize=1024;
  char buffer[bsize];
  
  size_t nm=3;
  regmatch_t match[nm];
  char *bra; // current position in read buffer
  char *ket; // position after a number was read successfully
  int end_block_match;
  value_t value;
  var_ndim_array *nda;

  nda=nda_alloc();
  if (target_type==DOUBLE_BLOCK)
    nda->is_double=1;
  else if ((target_type==INTEGER_BLOCK))
    nda->is_int=1;
  else
    exit(-5);

  do{ // read block line by line, each line is read by repeated calls
      // to strtod or strtol. lines and columsn are counted while reading
    do{ // get lines until it is not blank (or just a comment)
      fgets(buffer,bsize,f);
      //printf("%s (is comment: %i)\n",buffer,regexec(comment,buffer,0,NULL,0)==0);
    } while (regexec(&comment[0],buffer,0,NULL,0)==0 || strlen(buffer)<2);
    // check whether the line closes with a comment
    // if so, terminate the string at beginning of comment ⟨ ⟩
    if (regexec(&comment[1],buffer,nm,match,0)==0) buffer[match[1].rm_so]='\0';
    end_block_match=regexec(end,buffer,0,NULL,0);
    if (end_block_match!=0) {
      if (r>1 && c!=pre_c){ // beginning with the second row, we check
			    // the number of columns for consistency
	fprintf(stderr,"[read_block_nd] inconsistent number of columns, row=%i\n",r);
	exit(-3);
      } else {
	pre_c=c;
	c=0;
      }
      ket=buffer;
      printf("#⟨");
      if (target_type==DOUBLE_BLOCK){
	do {
	  bra=ket;
	  value.d=strtod(bra,&ket);
	  if (bra!=ket){
	    printf(" %g ",value.d);
	    c++;
	    nda_push(nda,&value);
	  }
	} while (bra!=ket);
      } else if (target_type==INTEGER_BLOCK){
	do {
	  bra=ket;
	  value.i=strtol(bra,&ket,10);
	  if (bra!=ket){
	    printf(" %i ",value.i);
	    c++;
	    nda_push(nda,&value);
	  }
	} while (bra!=ket);
      } else {
	fprintf(stderr,"[read_block_nd] error for target type: %i\n",target_type);
	exit(-3);
      }
      if (c>0) {
	r++; // so, this row contained at least one number.
      }
      printf("⟩\n");
    }
  } while (end_block_match!=0);
  printf("#");  printf("r=%i, c=%i, ",r,c);
  //adjust the dimension information in the n-dimensional array
  nda->nd=(r>1)?2:1;
  printf("nd=%i.\n",nda->nd);
  nda->size=(int *) malloc(sizeof(int)*(nda->nd));
  if(r>1){
    nda->size[0]=r;
    nda->size[1]=c;
  } else {
    nda->size[0]=c;
  }
  if (nda->N != r*c){
    fprintf(stderr,"[read_block_nd] inconsistent number of data-points: %i!=(%i×%i)\n",nda->N,r,c);
    exit(-3);
  }
  fflush(stdout);
  return nda;
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




/* reads data blocks from cfg file. The data is stored in dynamic
 * arrays first, then is copied into gsl data structures and is then
 * copied from there into a more convenient format.
 */

int read_problem_definition(FILE *cnf, ode_model_parameters *omp, const field_expression *fe, const regex_t *comment, main_options *cnf_options){
  int i,r,c;
  int C=0,T=0;
  int bsize=2048;
  char buffer[bsize];
  char *var_name;
  char *var_value;
  char *var_value_newline;
  var_ndim_array *nda;
  gsl_object *G;
  gsl_matrix_sd Data;
  gsl_matrix_sd ReferenceData;
 
  //char *c;
  //  size_t nm=3;
  //  regmatch_t match[nm];
  omp->normalisation_type=DATA_IS_ABSOLUTE;
  const field_expression *current;


  G=gsl_object_alloc();
  do {
    fgets(buffer,bsize,cnf);    
    // discard all comment lines
    if (regexec(&comment[0],buffer,0,NULL,0)!=0 ){
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
	
	if (strcmp(cnf_options->output_file,"sample.dat")==0 && strcmp(buffer,"output")==0)
	  strcpy(cnf_options->output_file,var_value);
	else if (cnf_options->sample_size<0 && strcmp(buffer,"sample_size")==0)
	  cnf_options->sample_size=strtol(var_value,NULL,0);
	else if (cnf_options->target_acceptance<0 && strcmp(buffer,"acceptance")==0)
	  cnf_options->target_acceptance=strtod(var_value,NULL);
	else if (cnf_options->initial_stepsize<0 && strcmp(buffer,"step_size")==0)
	  cnf_options->initial_stepsize=strtod(var_value,NULL);
	else if (strcmp(buffer,"C")==0){
	  C=strtod(var_value,NULL);
	  omp->size->C=C;
	  omp->E=(experiment**) malloc(sizeof(experiment*)*C);
	  for (i=0;i<C;i++) {
	    omp->E[i]=(experiment*) malloc(sizeof(experiment));
	  }
	  omp->ref_E=(experiment*) malloc(sizeof(experiment)); // just in case
	  printf("# %i experiment structures allocated.\n",C);
	}
	else if (strcmp(buffer,"t0")==0) {
	  omp->t0=strtod(var_value,NULL);
	  printf("# t0 = %g\n",omp->t0);
	}
	fflush(stdout);
      } else if (strchr(buffer,'[')!=NULL && strchr(buffer,']')!=NULL) {
	// so, this could be a tag.
	current=fe; // reset field expression (dynamic array) pointer to top
	while (current!=NULL && regexec(current->opening_bracket,buffer,0,NULL,0)!=0){
	  //printf("id_%i ",current->id);	
	  current=current->next;
	}
	printf("\n");
	printf("# found match with regular expression id %i: %s; ",current->id,buffer);
	switch (current->id){// reading the data
	case i_time:
	  nda=read_block_nd(cnf,DOUBLE_BLOCK,current->closing_bracket,comment);
	  gsl_object_init(G);
	  nda_to_gsl(nda,G);
	  r=G->gsl->matrix->size1;
	  c=G->gsl->matrix->size2;
	  for (i=0;i<C;i++){
	    omp->E[i]->t=gsl_vector_alloc(c);
	    gsl_matrix_get_row(omp->E[i]->t, G->gsl->matrix, i%r);
	  }
	  nda_free(nda);
	  printf("# measurement time(s) read.\n");	  
	  break;
	case i_reference_input:
	  //printf("# reading reference input.\n");
	  nda=read_block_nd(cnf,DOUBLE_BLOCK,current->closing_bracket,comment);
	  gsl_object_init(G);
	  nda_to_gsl(nda,G);
	  r=G->gsl->matrix->size1;
	  c=G->gsl->matrix->size2;
	  omp->ref_E->input_u=gsl_vector_alloc(c);
	  gsl_matrix_get_row(omp->ref_E->input_u,G->gsl->matrix,0);
	  nda_free(nda);
	  //read_block(1,ps->U,cnf,DOUBLE_BLOCK,omp->reference_u->data,comment);
	  printf("# reference input read.\n");
	  break;
	case i_reference_data: // reference data is used in
	  // normalisation, so we set it aside
	  // for now.
	  omp->normalisation_type=DATA_NORMALISED_BY_REFERENCE;
	  nda=read_block_nd(cnf,DOUBLE_BLOCK,current->closing_bracket,comment);
	  gsl_object_init(G);
	  nda_to_gsl(nda,G);
	  r=G->gsl->matrix->size1;
	  c=G->gsl->matrix->size2;
	  ReferenceData.M=gsl_matrix_alloc(r,c);
	  gsl_matrix_memcpy(ReferenceData.M,G->gsl->matrix);
	  omp->size->T=r; // number of time points for normalised_by_reference data
	  nda_free(nda);
	  printf("# reference data read.\n");
	  break;
	case i_sd_reference_data:
	  nda=read_block_nd(cnf,DOUBLE_BLOCK,current->closing_bracket,comment);
	  gsl_object_init(G);
	  nda_to_gsl(nda,G);
	  r=G->gsl->matrix->size1;
	  c=G->gsl->matrix->size2;
	  ReferenceData.sd=gsl_matrix_alloc(r,c);
	  gsl_matrix_memcpy(ReferenceData.sd, G->gsl->matrix);
	  nda_free(nda);	  
	  //read_block(ps->T,ps->F,cnf,DOUBLE_BLOCK,RD->sd->data,comment);
	  printf("# standard deviation of reference data read.\n");
	  break;
	case i_input:
	  nda=read_block_nd(cnf,DOUBLE_BLOCK,current->closing_bracket,comment);
	  gsl_object_init(G);
	  nda_to_gsl(nda,G);
	  r=G->gsl->matrix->size1;
	  c=G->gsl->matrix->size2;
	  omp->size->U=c;
	  for (i=0;i<C;i++){
	    omp->E[i]->input_u=gsl_vector_alloc(c);
	    gsl_matrix_get_row(omp->E[i]->input_u,G->gsl->matrix,i%r);
	  }
	  nda_free(nda);	  
	  //read_block(ps->C,ps->U,cnf,DOUBLE_BLOCK,omp->input_u->data,comment);
	  printf("# input read.\n");
	  break;
	case i_initial_conditions:
	  nda=read_block_nd(cnf,DOUBLE_BLOCK,current->closing_bracket,comment);
	  gsl_object_init(G);
	  nda_to_gsl(nda,G);
	  r=G->gsl->matrix->size1;
	  c=G->gsl->matrix->size2;
	  for (i=0;i<C;i++){
	    omp->E[i]->init_y=gsl_vector_alloc(c);
	    gsl_matrix_get_row(omp->E[i]->init_y,G->gsl->matrix,i%r);
	  }
	  nda_free(nda);	  
	  //	  read_block(ps->C,ps->N,cnf,DOUBLE_BLOCK,omp->initial_conditions_y->data,comment);
	  printf("# initial conditions read.\n");
	  break;
	case i_data: // data possibly needs to be normalised, so we
	  // set it aside for now;
	  nda=read_block_nd(cnf,DOUBLE_BLOCK,current->closing_bracket,comment);
	  gsl_object_init(G);
	  nda_to_gsl(nda,G);
	  r=G->gsl->matrix->size1;
	  c=G->gsl->matrix->size2;
	  Data.M=gsl_matrix_alloc(r,c);
	  gsl_matrix_memcpy(Data.M,G->gsl->matrix);
	  nda_free(nda);
	  //read_block((ps->C)*(ps->T),ps->F,cnf,DOUBLE_BLOCK,omp->Data->data,comment);
	  printf("# data read.\n");
	  break;
	case i_sd_data:
	  nda=read_block_nd(cnf,DOUBLE_BLOCK,current->closing_bracket,comment);
	  gsl_object_init(G);
	  nda_to_gsl(nda,G);
	  r=G->gsl->matrix->size1;
	  c=G->gsl->matrix->size2;
	  Data.sd=gsl_matrix_alloc(r,c);
	  gsl_matrix_memcpy(Data.sd, G->gsl->matrix);
	  nda_free(nda);
	  //read_block((ps->C)*(ps->T),ps->F,cnf,DOUBLE_BLOCK,omp->sdData->data,comment);
	  printf("# standard deviation of data read.\n");
	  break;
	case i_prior_mu:
	  nda=read_block_nd(cnf,DOUBLE_BLOCK,current->closing_bracket,comment);
	  gsl_object_init(G);
	  nda_to_gsl(nda,G);
	  r=G->gsl->matrix->size1;
	  c=G->gsl->matrix->size2;
	  omp->size->D=(r>c?r:c);
	  omp->prior_mu=gsl_vector_alloc(omp->size->D);
	  memcpy(omp->prior_mu->data,G->gsl->matrix->data,omp->size->D*sizeof(double));
	  nda_free(nda);
	  //read_block(ps->D,1,cnf,DOUBLE_BLOCK,omp->prior_mu->data,comment);
	  printf("# prior mean read.\n");
	  break;
	case i_prior_icov:
	  nda=read_block_nd(cnf,DOUBLE_BLOCK,current->closing_bracket,comment);
	  gsl_object_init(G);
	  nda_to_gsl(nda,G);
	  r=G->gsl->matrix->size1;
	  c=G->gsl->matrix->size2;
	  omp->prior_inverse_cov=gsl_matrix_alloc(r,c);
	  gsl_matrix_memcpy(omp->prior_inverse_cov, G->gsl->matrix);
	  nda_free(nda);
	  //read_block(ps->D,ps->D,cnf,DOUBLE_BLOCK,omp->prior_inverse_cov->data,comment);
	  printf("# prior inverse covariance matrix read.\n");
	  break;
	case i_norm_f: // norm_f and norm_t can be combined, so they add up.
	  omp->normalisation_type++;
	  nda=read_block_nd(cnf,DOUBLE_BLOCK,current->closing_bracket,comment);
	  gsl_object_init(G);
	  nda_to_gsl(nda,G);
	  r=G->gsl->matrix_int->size1;
	  c=G->gsl->matrix_int->size2;
	  omp->norm_f=gsl_matrix_int_alloc(r,c);
	  gsl_matrix_int_memcpy(omp->norm_f, G->gsl->matrix_int);
	  nda_free(nda);
	  //read_block(ps->C,ps->F,cnf,INTEGER_BLOCK,omp->norm_f->data,comment);
	  printf("# norm_f read.\n");
	  break;
	case i_norm_t: // mixing reference data and these two types will then lead to an error.
	  omp->normalisation_type+=2;
	  nda=read_block_nd(cnf,DOUBLE_BLOCK,current->closing_bracket,comment);
	  gsl_object_init(G);
	  nda_to_gsl(nda,G);
	  r=G->gsl->matrix_int->size1;
	  c=G->gsl->matrix_int->size2;
	  omp->norm_t=gsl_matrix_int_alloc(r,c);
	  gsl_matrix_int_memcpy(omp->norm_t, G->gsl->matrix_int);
	  nda_free(nda);
	  //	  read_block(ps->C,ps->F,cnf,INTEGER_BLOCK,omp->norm_t->data,comment);
	  printf("# norm_t read.\n");
	  break;
	default:
	  fprintf(stderr,"unknown tag: %s\n",buffer);
	  exit(-1);
	  break;
	}// switch what to do.
	//printf("/switch\n");
      }// if var_name=var_value
      //printf("/if\n");
    }
  } while (!feof(cnf));
  fflush(stderr);
  r=0;
  for (i=0;i<C;i++){
    T=omp->E[i]->t->size;
    r+=T;
  }
  //printf("cumulative T: CT=%i\n",r);fflush(stdout);
  omp->Data=gsl_matrix_alloc(r,omp->size->F);
  omp->sdData=gsl_matrix_alloc(r,omp->size->F);
  ode_model_parameters_alloc(omp);
  ode_model_parameters_link(omp); // all "views" are linked to their data objects
  fflush(stdout);
  if (omp->normalisation_type>0){
    printf("# calculating relative data.\n");
    printf("# normalisation type: %i.\n",omp->normalisation_type);
    switch (omp->normalisation_type){
    case DATA_NORMALISED_BY_REFERENCE:
      printf("# taking ratios of Data and Refrence Data...");
      // all experiments, must be congruent
      // i.e.: omp->E[i]->t must all match, as well as all F
      // this can be done, by simply having just one time row
      ratio_with_sd(&Data,&ReferenceData); //Data/ReferenceData.
      gsl_matrix_memcpy(omp->Data,Data.M);
      gsl_matrix_memcpy(omp->sdData,Data.sd);
      printf("# done.\n");
      break;
    case DATA_NORMALISED_BY_TIMEPOINT:
      printf("# normalising using one time instance... %i\n",gsl_matrix_int_get(omp->norm_t,0,0));
      fflush(stdout);
      gsl_matrix_memcpy(omp->Data,Data.M);
      gsl_matrix_memcpy(omp->sdData,Data.sd);
      normalise_by_timepoint_with_sd(omp);
      printf("# done.\n");
      break;
    case DATA_NORMALISED_BY_STATE_VAR:
      printf("# normalising using another state variable at specified time instance...");
      gsl_matrix_memcpy(omp->Data,Data.M);
      gsl_matrix_memcpy(omp->sdData,Data.sd);
      normalise_by_state_var_with_sd(omp);
      printf("# done.\n");
      break;
    }
  }
  // cleanup
  gsl_matrix_free(Data.M);
  gsl_matrix_free(Data.sd);
  gsl_matrix_free(ReferenceData.M);
  gsl_matrix_free(ReferenceData.sd);
  fflush(stdout);
  gsl_object_free(G);
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
  strcpy(fn->name[i_initial_conditions],"init(ial_conditions)?");
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

int parse_config(FILE *cnf, ode_model_parameters *omp,  main_options *cnf_options){
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
  fe=field_expression_init(&fn);
  /*comments in file: a line that is only a comment: "_blank_ #
   * comment" must be disregarded when line counting; a line that
   * contains both content and comments must be read up to the comment
   * symbol.
   */
  regcomp(&comment[0],"^\\s*((#|//|%)|$)",REG_EXTENDED); // a blank line or a comment line
  regcomp(&comment[1],"\\w+\\s*(#|//|%)",REG_EXTENDED); // some content and a comment
  
  omp->t0=cnf_options->t0;
  read_problem_definition(cnf, omp, fe, comment, cnf_options);
  printf("# configuration read.\n");
  fflush(stdout);
  while (fe!=NULL) { //empty stack
    regfree(fe->opening_bracket);
    regfree(fe->closing_bracket);
    fe=fe->next;    
  }
  fflush(stdout);

  return EXIT_SUCCESS;
}
