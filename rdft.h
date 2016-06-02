#ifndef MYRDFT_H
#define MYRDFT_H

#include "test.h"

void rdft_original_slow(double *a, double *b, double *x, double *xi, double *xia);
void dft_matrix_complex_double(dcomplex *f);
void r_matrix_complex_double(dcomplex *r);

void perm_matrix_double(double *r);

#endif
