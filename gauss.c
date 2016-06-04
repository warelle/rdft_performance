#include "matrix_double.h"
#include "gauss.h"
#include "gen.h"
#include "lu.h"
#include <string.h>

#include <stdio.h>

void gauss_matrix(double *g){
  for(int i=0; i<MATRIX_SIZE; i++)
    for(int j=0; j<MATRIX_SIZE; j++)
      IDX(g, MATRIX_SIZE, i, j) = rand_normal_double(0,1);
}

void solve_with_gauss_iteration_double(double *a, double *b, double *x, double *xi, double *xia){
  double *g, *ag, *y;
  double alpha=1.0, zero=0.0;
  int size = MATRIX_SIZE, inc=1;
  char non = 'N';
  char l='L', u='U';

  alloc_matrix_double(&g, MATRIX_SIZE);
  alloc_matrix_double(&ag, MATRIX_SIZE);
  alloc_vector_double(&y, MATRIX_SIZE);
  gauss_matrix(g);

  // ag <= A*G
  dgemm_(&non,&non, &size,&size,&size, &alpha, a,&size, g,&size, &zero, ag, &size);

  // lu steps
  dgetrfw_(&size, &size, ag, &size);

  // back-forward
  dcopy_(&size, b,&inc, y,&inc);

  dtrsm_(&l,&l,&non, &u,   &size,&inc, &alpha, ag, &size, y, &size);
  dtrsm_(&l,&u,&non, &non, &size,&inc, &alpha, ag, &size, y, &size);

  // x = Gx
  dgemv_(&non, &size, &size, &alpha, g, &size, y, &inc, &zero, x, &inc);

  // iteration_double(d_fa, d_l, d_u, d_fb, xi);
  // iteration_double_another(d_f, a, d_l, d_u, b, xia);

  free_matrix_double(&g);
  free_matrix_double(&ag);
  free_vector_double(&y);
}
