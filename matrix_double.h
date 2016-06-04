#ifndef MYDMAT_H
#define MYDMAT_H

#include "global.h"

// ----- alloc/free ----- //
#include <stdlib.h>

// double
static void alloc_vector_double(double **a, int size){
	(*a) = (double*)malloc(sizeof(double)*size);
}
static void free_vector_double(double **a){
	free(*a);
}
static void alloc_matrix_double(double **a, int size){
	(*a) = (double*)malloc(sizeof(double)*size*size);
}
static void free_matrix_double(double **a){
	free(*a);
}

#endif
