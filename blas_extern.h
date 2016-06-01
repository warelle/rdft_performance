#ifndef MYBLAS_EXTERN_H
#define MYBLAS_EXTERN_H

//----------------- LEVEL 1 BLAS ------------------//

// dim, X[dim], x_inc, Y[dim], y_inc
// Y = X
extern void dcopy_(int*, double*,int*, double*,int*);

// dim, alpha, X[dim], x_inc, Y[dim], y_inc
//  Y = Y + alpha*X
extern void daxpy_(int*, double*, double*, int*, double*, int*);

// dim, X[dim], x_inc
// return norm(X)
extern double dnrm2_(int *n, double *x, int *incx);

//----------------- LEVEL 2 BLAS ------------------//
extern void dgemv_(char *trans, int *m, int *n, double *alpha, double *A, int *ldA, double *x, int *incx, double *beta , double *y, int *incy);

//----------------- LEVEL 3 BLAS ------------------//
extern void dgemm_(char*, char*, int*, int*,int*, double*, double*, int*, double*, int*, double*, double*, int*);


//----------------- LAPACK ------------------//
extern void dsgesv_(int*,int*,double*,int*,int*,double*,int*,double*,int*,double*,float*,int*,int*);
extern void dtrsm_ (char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);


#endif
