#ifndef MYGALLOC_H
#define MYGALLOC_H

#include "global.h"
#include "givens.h"

// ----- alloc/free ----- //
#include <stdlib.h>
#include <stdio.h>

// matrix int
static void alloc_matrix_int(int **a, int size){
	(*a) = (int*)malloc(sizeof(int)*size*size);
  if(*a == NULL){
    fprintf(stderr,"malloc failed\n");
    exit(0);
  }
}
static void free_matrix_int(int **a){
	free(*a);
  *a = NULL;
}

#endif
