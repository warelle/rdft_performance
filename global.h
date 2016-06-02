#ifndef MYGLOBAL_H
#define MYGLOBAL_H

// ----- C99 complex ----- //
#include <complex.h>
typedef double _Complex dcomplex;
#define CNUM(a,b) ((a)+I*(b))

// ----- matrix access ----- //
#define IDX(mat,n, i,j) (mat)[(i)+(j)*(n)]

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
// dcomplex
static void alloc_vector_complex_double(dcomplex **a, int size){
	(*a) = (dcomplex*)malloc(sizeof(dcomplex)*size);
}
static void free_vector_complex_double(dcomplex **a){
	free(*a);
}
static void alloc_matrix_complex_double(dcomplex **a, int size){
	(*a) = (dcomplex*)malloc(sizeof(dcomplex)*size*size);
}
static void free_matrix_complex_double(dcomplex **a){
	free(*a);
}

#endif
