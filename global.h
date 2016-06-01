#ifndef MYGLOBAL_H
#define MYGLOBAL_H

// ----- C99 complex ----- //
#include <complex.h>
typedef double _Complex dcomplex;

// ----- matrix access ----- //
#define IDX(mat,n, i,j) (mat)[(i)*(n)+(j)]

// ----- alloc/free ----- //
#include <stdlib.h>
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
