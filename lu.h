#ifndef MYLU_H
#define MYLU_H

#include <math.h>
#include <complex.h>
#include "test.h"

// LU decomposition
//  input : m
void lu_double(double *m);
//void lu_complex_double(std::complex<double> m[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE]);

// LU solver
//void l_step_complex_double(std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> y[MATRIX_SIZE]);
//void u_step_complex_double(std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> y[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE]);


#endif
