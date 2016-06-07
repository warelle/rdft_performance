#ifndef MYBLAS_EXTERN_H
#define MYBLAS_EXTERN_H

#include "global.h"

//----------------- LEVEL 1 BLAS ------------------//

// dim, X[dim], x_inc, Y[dim], y_inc
// Y = X
extern void dcopy_(int*, double*,int*, double*,int*);
extern void zcopy_(int*, dcomplex*,int*, dcomplex*,int*);

// dim, alpha, X[dim], x_inc, Y[dim], y_inc
//  Y = Y + alpha*X
extern void daxpy_(int*, double*, double*, int*, double*, int*);
extern void zaxpy_(int*, dcomplex*, dcomplex*, int*, dcomplex*, int*);

// dim, X[dim], x_inc
// return norm(X)
extern double dnrm2_(int*, double*, int*);
extern double znrm2_(int*, dcomplex*, int*);

extern void drot_(int*, double*, int*, double*, int*, double*, double*);

extern void zscal_(int*, dcomplex*, dcomplex*, int*);

//----------------- LEVEL 2 BLAS ------------------//
extern void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double* , double*, int*);
extern void zgemv_(char*, int*, int*, dcomplex*, dcomplex*, int*, dcomplex*, int*, dcomplex* , dcomplex*, int*);

//----------------- LEVEL 3 BLAS ------------------//
extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

extern void zgemm_(char*, char*, int*, int*, int*, dcomplex*, dcomplex*, int*, dcomplex*, int*, dcomplex*, dcomplex*, int*);

//----------------- LAPACK ------------------//
extern void dsgesv_(int*,int*,double*,int*,int*,double*,int*,double*,int*,double*,float*,int*,int*);
extern void dtrsm_ (char*, char*, char*, char*, int*, int*, double*, double*, int*, double*, int*);

extern void zsgesv_(int*,int*,dcomplex*,int*,int*,dcomplex*,int*,dcomplex*,int*,dcomplex*,float*,int*,int*);
extern void ztrsm_ (char*, char*, char*, char*, int*, int*, dcomplex*, dcomplex*, int*, dcomplex*, int*);


//----------------- LAPACK/mod ------------------//
extern void dgetrfw_(int*, int*, double*, int*);
extern void zgetrfw_(int*, int*, dcomplex*, int*);


#endif
