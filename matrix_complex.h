#ifndef MYCMAT_H
#define MYCMAT_H

#include "global.h"

// ----- alloc/free ----- //
#include <stdlib.h>

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
