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

void solve_with_gauss_iteration_double(double *a, double *b, double *x, double *xi, double *xia, double *g, double *ag){
  double alpha=1.0, zero=0.0;
  int size = MATRIX_SIZE, inc=1;

  gauss_matrix(g);

  // ag <= A*G
  dgemm_("N","N", &size,&size,&size, &alpha, a,&size, g,&size, &zero, ag, &size);

  // lu steps
  dgetrfw_(&size, &size, ag, &size);

  // back-forward
  dtrsm_("L","L","N","U", &size,&inc, &alpha, ag, &size, b, &size);
  dtrsm_("L","U","N","N", &size,&inc, &alpha, ag, &size, b, &size);

  // x = Gx
  dgemv_("N", &size, &size, &alpha, g, &size, b, &inc, &zero, x, &inc);

  // iteration_double(d_fa, d_l, d_u, d_fb, xi);
  // iteration_double_another(d_f, a, d_l, d_u, b, xia);
}
