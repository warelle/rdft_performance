#ifndef MYCMAT_H
#define MYCMAT_H

#include "global.h"

// ----- alloc/free ----- //
#include <stdlib.h>
#include <stdio.h>

// dcomplex
static void alloc_vector_complex_double(dcomplex **a, int size){
	(*a) = (dcomplex*)malloc(sizeof(dcomplex)*size);
  if(*a == NULL){
    fprintf(stderr,"malloc failed\n");
    exit(0);
  }
}
static void free_vector_complex_double(dcomplex **a){
	free(*a);
  *a = NULL;
}
static void alloc_matrix_complex_double(dcomplex **a, int size){
	(*a) = (dcomplex*)malloc(sizeof(dcomplex)*size*size);
  if(*a == NULL){
    fprintf(stderr,"malloc failed\n");
    exit(0);
  }
}
static void free_matrix_complex_double(dcomplex **a){
	free(*a);
  *a = NULL;
}

#endif
