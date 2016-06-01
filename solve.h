#ifndef MYSOLVE_H
#define MYSOLVE_H

#include "test.h"

// solve Ax = b
// input : a,b
// return: x

// solver
//void solve_with_rdft_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE]);
//void solve_with_rdft_perm_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE], int *counter=NULL);
//void solve_with_rdft_givens_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE], int *counter=NULL);
//void solve_with_rdft_givens_two_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE], int *counter=NULL);
//void solve_with_rdft_both_givens_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE], int *counter=NULL);
//
//void solve_without_pivot_double(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double x[MATRIX_SIZE], double xi[MATRIX_SIZE], double xia[MATRIX_SIZE]);
void solve_with_partial_pivot(double *a, double *b, double *x, double *xi, double *xia);

// iteration
// void iteration_double(double a[MATRIX_SIZE][MATRIX_SIZE], double l[MATRIX_SIZE][MATRIX_SIZE], double u[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double x[MATRIX_SIZE]);
// void iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE]);
//
// void iteration_double_another(double f[MATRIX_SIZE][MATRIX_SIZE], double a[MATRIX_SIZE][MATRIX_SIZE], double l[MATRIX_SIZE][MATRIX_SIZE], double u[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double x[MATRIX_SIZE]);
// void iteration_complex_double_another(std::complex<double> f[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE]);

#endif
