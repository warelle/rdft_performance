#include "lu.h"

#include <stdio.h>

void lu_double(double *m){
  int i,j,k;

  for(i=0; i<MATRIX_SIZE-1; i++){
    for(j=i+1; j<MATRIX_SIZE; j++){
      IDX(m,MATRIX_SIZE,j,i) /= IDX(m,MATRIX_SIZE,i,i);
      for(k=i+1; k<MATRIX_SIZE; k++){
        IDX(m,MATRIX_SIZE,j,k) -= IDX(m,MATRIX_SIZE,j,i)*IDX(m,MATRIX_SIZE,i,k);
      }
    }
  }
}
/*
void lu_complex_double(dcomplex *m){
  int i,j,k;

  for(i=0; i<MATRIX_SIZE-1; i++){
    for(j=i+1; j<MATRIX_SIZE; j++){
      mat[j][i] /= mat[i][i];
      for(k=i+1; k<MATRIX_SIZE; k++){
        mat[j][k] -= mat[j][i]*mat[i][k];
      }
    }
  }
}
*/
/*
void l_step_complex_double(std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> y[MATRIX_SIZE]){
  int i,j;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
      y[i] = b[i];

  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<i; j++){
      y[i] -= l[i][j]*y[j];
    }
    y[i] /= l[i][i];
  }
}
void u_step_complex_double(std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> y[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE]){
  int i,j;

  // [TODO] more efficient code
  for(i=0; i<MATRIX_SIZE; i++)
      x[i] = y[i];

  for(i=MATRIX_SIZE-1; i>=0; i--){
    for(j=MATRIX_SIZE-1; j>=i+1; j--){
      x[i] -= u[i][j]*x[j];
    }
    x[i] /= u[i][i];
  }
}
*/
