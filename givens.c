#include "givens.h"
#include "gen.h"
#include "givens_alloc.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef THREAD_NUM
#define GS_NUM THREAD_NUM
#else
#define GS_NUM 1
#endif

#define MATH_PI 3.14159265358979

// ----- create/delete givens matrix basic ----- //
// givens matrix
givens_matrix *create_givens_matrix(int i, int j, double theta){
  givens_matrix *gm = (givens_matrix*)malloc(sizeof(givens_matrix));
  if(gm == NULL){
    fprintf(stderr,"malloc failed\n");
    exit(0);
  }

  if(i > j){
    int tmp = i;
    i = j;
    j = tmp;
  }

  gm->i = i;
  gm->j = j;
  gm->c = cos(theta);
  gm->s = sin(theta);

  return gm;
}
void delete_givens_matrix(givens_matrix *gm){
  if(gm) free(gm);
  gm = NULL;
}

// givens matrix list
givens_matrix_list *create_givens_matrix_list(){
  givens_matrix_list *gml = (givens_matrix_list*)malloc(sizeof(givens_matrix_list));
  if(gml == NULL){
    fprintf(stderr,"malloc failed\n");
    exit(0);
  }
  gml->gm = NULL;
  gml->next = NULL;
  return gml;
}
void delete_givens_matrix_list(givens_matrix_list *gml){
  if(gml){
    if(gml->next != NULL){
      delete_givens_matrix_list(gml->next);
    }
    delete_givens_matrix(gml->gm);
    free(gml);
    gml = NULL;
  }
}

// givens sequence list
givens_sequence_list *create_givens_sequence_list(){
  givens_sequence_list *gsl,*cur;
  for(int i=0; i<GS_NUM; i++){
    givens_sequence_list *tmp = (givens_sequence_list*)malloc(sizeof(givens_sequence_list));
    if(tmp == NULL){
      fprintf(stderr,"malloc failed\n");
      exit(0);
    }
    tmp->gml  = NULL;
    tmp->next = NULL;
    cur->next = tmp;
    cur = tmp;

    // assign top
    if(i==0) gsl = tmp;
  }
  return gsl;
}
void delete_givens_sequence_list(givens_sequence_list *gsl){
  if(gsl){
    if(gsl->next != NULL){
      delete_givens_sequence_list(gsl->next);
    }
    free(gsl);
    gsl = NULL;
  }
}

// ----- create/delete givens matrix ----- //
givens_sequence_list *generate_givens(){
  givens_sequence_list *gsl_r, *gsl_cur;
  givens_matrix_list *gml=NULL, *gml_cur;
  int *shfl;

  alloc_matrix_int(&shfl, MATRIX_SIZE);

  for(int i=0; i<MATRIX_SIZE; i++)
    shfl[i] = i;

  // Fisher-Yates shuffle
  for(int i=0; i<MATRIX_SIZE-1; i++){
    int k = uniform_int(i,MATRIX_SIZE-1);
    int tmp = shfl[k];
    shfl[k] = shfl[i];
    shfl[i] = tmp;
  }

  gsl_r = create_givens_sequence_list();
  gsl_cur = gsl_r;
  const int gm_num = MATRIX_SIZE/GS_NUM;
  for(int i=0; i<GS_NUM; i++){
    gsl_cur = create_givens_sequence_list();
    for(int j=0; j<gm_num/2; j++){
      gml_cur = create_givens_matrix_list();
      givens_matrix *gm = create_givens_matrix(shfl[gm_num*i + 2*j],shfl[gm_num*i+2*j+1], uniform()*2*MATH_PI);
      gml_cur->gm = gm;

      if(j == 0){
        gsl_cur->gml = gml_cur;
      }else{
        gml->next = gml_cur;
      }
      gml = gml_cur;
    }
    if(i == 0)
      gsl_r = gsl_cur;
    gsl_cur = gsl_cur->next;
  }

  free_matrix_int(&shfl);

  return gsl_r;
}
void delete_givens(givens_sequence_list *gsl){
  delete_givens_sequence_list(gsl);
}



// ----- multipliation methods ----- //
// ----- givens matrix-vector multiplication ----- //
void mat_vec_dot_givens(givens_matrix *gm, double *a){
  int inc = 1;
  gm->s = -gm->s;
  drot_(&inc, a+(gm->i),&inc, a+(gm->j),&inc, &(gm->c), &(gm->s));
  gm->s = -gm->s;
}
void mat_vec_dot_givens_matrix_list(givens_matrix_list *gml, double *a){
  while(gml != NULL){
    mat_vec_dot_givens(gml->gm, a);
    gml = gml->next;
  }
}
void mat_vec_dot_givens_sequence(givens_sequence_list *gsl, double *a){
  while(gsl != NULL){
    mat_vec_dot_givens_matrix_list(gsl->gml, a);
    gsl = gsl->next;
  }
}

// ----- givens matrix-matrix multiplication ----- //
void mat_mul_givens(givens_matrix *gm, double *a){
  int dim = MATRIX_SIZE, inc = 1;
  drot_(&dim, a+(MATRIX_SIZE*(gm->i)),&inc, a+(MATRIX_SIZE*(gm->j)),&inc, &(gm->c), &(gm->s));
}
void mat_mul_givens_matrix_list(givens_matrix_list *gml, double *a){
  while(gml != NULL){
    mat_mul_givens(gml->gm, a);
    gml = gml->next;
  }
}
void mat_mul_givens_sequence(givens_sequence_list *gsl, double *a){
  while(gsl != NULL){
    mat_mul_givens_matrix_list(gsl->gml, a);
    gsl = gsl->next;
  }
}
