#ifndef MYDMAT_H
#define MYDMAT_H

#include "global.h"

// ----- alloc/free ----- //
#include <stdlib.h>
#include <stdio.h>

// double
static void alloc_vector_double(double **a, int size){
	(*a) = (double*)malloc(sizeof(double)*size);
  if(*a == NULL){
    fprintf(stderr,"malloc failed\n");
    exit(0);
  }
}
static void free_vector_double(double **a){
	free(*a);
  *a = NULL;
}
static void alloc_matrix_double(double **a, int size){
	(*a) = (double*)malloc(sizeof(double)*size*size);
  if(*a == NULL){
    fprintf(stderr,"malloc failed\n");
    exit(0);
  }
}
static void free_matrix_double(double **a){
	free(*a);
  *a = NULL;
}

#endif
