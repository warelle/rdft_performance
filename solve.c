#include "solve.h"
#include "lu.h"
#include <stdio.h>
#include <stdlib.h>

// solving
double *d_f;
double *d_r;
double *d_fr;
double *d_fra;
double *d_l;
double *d_u;
double *d_y;
double *d_frb;

double d_w[MATRIX_SIZE];
double d_v[MATRIX_SIZE];
double d_z[MATRIX_SIZE];
double d_ax[MATRIX_SIZE];

//std::complex<double> cd_w[MATRIX_SIZE];
//std::complex<double> cd_v[MATRIX_SIZE];
//std::complex<double> cd_z[MATRIX_SIZE];
//std::complex<double> cd_ax[MATRIX_SIZE];

void solve_no_pivoting(double *a, double *b, double *x, double *xi, double *xia){
  double alpha=1.0;
  int size = MATRIX_SIZE, inc=1;

  // lu steps
  dgetrfw_(&size, &size, a, &size);

  // back-forward
  dcopy_(&size, b,&inc, x,&inc);

  dtrsm_("L","L","N","U", &size,&inc, &alpha, a, &size, x, &size);
  dtrsm_("L","U","N","N", &size,&inc, &alpha, a, &size, x, &size);
}

void solve_with_partial_pivot(double *a, double *b, double *x, double *xi, double *xia){
  int dim = MATRIX_SIZE;
  int nrhs = 1, inc = 1;
  double alpha = 1.0;
  int info;

  int *piv = (int*)malloc(sizeof(int)*MATRIX_SIZE);
  if(piv == NULL)
    fprintf(stderr, "malloc error in solve.c\n");

  dgesv_(&dim, &nrhs, a,&dim, piv, b,&dim, &info);

  dcopy_(&dim, b,&inc, x,&inc);

  free(piv);

  // iteration_double(d_ap,d_l, d_u, d_bp, xi);
  // iteration_double_another(d_perm, a,d_l, d_u, b, xia);
}

//void iteration_double(double a[MATRIX_SIZE][MATRIX_SIZE], double l[MATRIX_SIZE][MATRIX_SIZE], double u[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double x[MATRIX_SIZE]){
//  int itera = 0;
//  int iteramax = ITER_MAX;
//
//  while(itera++ < iteramax){
//    mat_vec_dot_double(a,x,d_ax);
//    vec_sub_double(b,d_ax,d_w);
//    l_step_double(l,d_w,d_v);
//    u_step_double(u,d_v,d_z);
//    vec_add_double(x,d_z,x);
//  }
//}
/*
void iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE]){
  int itera = 0;
  int iteramax = ITER_MAX;

  while(itera++ < iteramax){
    mat_vec_dot_complex_double(a,x,cd_ax);
    vec_sub_complex_double(b,cd_ax,cd_w);
    l_step_complex_double(l,cd_w,cd_v);
    u_step_complex_double(u,cd_v,cd_z);
    vec_add_complex_double(x,cd_z,x);
  }
}
*/
//void iteration_double_another(double f[MATRIX_SIZE][MATRIX_SIZE], double a[MATRIX_SIZE][MATRIX_SIZE], double l[MATRIX_SIZE][MATRIX_SIZE], double u[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE], double x[MATRIX_SIZE]){
//  int itera = 0;
//  int iteramax = ITER_MAX;
//
//  while(itera++ < iteramax){
//    mat_vec_dot_double(a,x,d_ax);
//    vec_sub_double(b,d_ax,d_w);
//    mat_vec_dot_double(f,d_w,d_ax);
//    l_step_double(l,d_ax,d_v);
//    u_step_double(u,d_v,d_z);
//    vec_add_double(x,d_z,x);
//  }
//}
/*
void iteration_complex_double_another(std::complex<double> f[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE]){
  int itera = 0;
  int iteramax = ITER_MAX;

  while(itera++ < iteramax){
    mat_vec_dot_complex_double(a,x,cd_ax);
    vec_sub_complex_double(b,cd_ax,cd_w);
    mat_vec_dot_complex_double(f,cd_w,cd_ax);
    l_step_complex_double(l,cd_ax,cd_v);
    u_step_complex_double(u,cd_v,cd_z);
    vec_add_complex_double(x,cd_z,x);
  }
}
*/



