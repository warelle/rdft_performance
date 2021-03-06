#ifndef MYRDFT_H
#define MYRDFT_H

#include "test.h"


void dft_matrix_complex_double(dcomplex *f);
void r_matrix_complex_double(dcomplex *r);

void rdft_original_slow(double *a, double *b, double *x, double *xi, double *xia);
void fftw_rdft_original(double *a, double *b, double *x, double *xi, double *xia, dcomplex *fra, dcomplex *r, dcomplex *frb);

void fftw_rdft_right_two_givens(double *a, double *b, double *x, double *xi, double *xia, dcomplex *fra, dcomplex *r, dcomplex *frb, double *ass);

#endif
