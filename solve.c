#include "solve.h"
#include "lu.h"
#include <stdio.h>

// solving
double *d_f;
double *d_r;
double *d_fr;
double *d_fra;
double *d_l;
double *d_u;
double *d_y;
double *d_frb;

//std::complex<double> cd_f[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_r[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_rp[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_rp_sub[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_rgivens[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_rgivens_sub[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_fr[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_frg[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_frg_tmp[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_fra[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_l[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_u[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_y[MATRIX_SIZE];
//std::complex<double> cd_frb[MATRIX_SIZE];
//std::complex<double> cd_tmp[MATRIX_SIZE][MATRIX_SIZE];
//
//std::complex<double> cd_r1[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_r2[MATRIX_SIZE][MATRIX_SIZE];
//std::complex<double> cd_r3[MATRIX_SIZE][MATRIX_SIZE];

double d_w[MATRIX_SIZE];
double d_v[MATRIX_SIZE];
double d_z[MATRIX_SIZE];
double d_ax[MATRIX_SIZE];

//std::complex<double> cd_w[MATRIX_SIZE];
//std::complex<double> cd_v[MATRIX_SIZE];
//std::complex<double> cd_z[MATRIX_SIZE];
//std::complex<double> cd_ax[MATRIX_SIZE];

void solve_with_partial_pivot(double *a, double *b, double *x, double *xi, double *xia){
  int dim = MATRIX_SIZE;
  int nrhs = 1;
  int iter,info;

  int *piv = (int*)malloc(sizeof(int)*MATRIX_SIZE);
  double *work  = (double*)malloc(sizeof(double)*MATRIX_SIZE);
  float  *workf = (float*) malloc(sizeof(float)*MATRIX_SIZE*MATRIX_SIZE);

  dsgesv_(&dim, &nrhs, a,&dim, piv, b,&dim, x,&dim, work, workf, &iter, &info);

  // iteration_double(d_ap,d_l, d_u, d_bp, xi);
  // iteration_double_another(d_perm, a,d_l, d_u, b, xia);
}

//void solve_without_pivot_double(double *a, double *b, double *x, double *xi, double *xia){
//  int dim = MATRIX_SIZE;
//  int nrhs = 1;
//
//
//  // iteration_double(d_ap,d_l, d_u, d_bp, xi);
//  // iteration_double_another(d_perm, a,d_l, d_u, b, xia);
//}

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

void solve_with_rdft_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE]){
  int i;

  dft_matrix_complex_double(cd_f);
  r_matrix_complex_double(cd_r);
  mat_mul_complex_double(cd_f,cd_r,cd_fr);
  mat_mul_complex_double(cd_fr,a,cd_fra);
  mat_vec_dot_complex_double(cd_fr,b,cd_frb);

  lu_complex_double(cd_fra, cd_l, cd_u);

  l_step_complex_double(cd_l, cd_frb, cd_y);
  u_step_complex_double(cd_u, cd_y, x);

  for(i=0; i<MATRIX_SIZE; i++){
    xi[i] = x[i];
    xia[i] = x[i];
  }

  iteration_complex_double(cd_fra, cd_l, cd_u, cd_frb, xi);
  iteration_complex_double_another(cd_fr, a, cd_l, cd_u, b, xia);
}

void solve_with_rdft_perm_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE], int *counter){
  int i;

  dft_matrix_complex_double(cd_f);
  r_perm_matrix_complex_double(cd_rp);
  mat_mul_complex_double(cd_f,cd_rp,cd_fr); // use cd_frg for tmp
  mat_mul_complex_double(cd_fr,a,cd_fra);
  mat_vec_dot_complex_double(cd_fr,b,cd_frb);

  if(counter != NULL){
    mat_mul_complex_double(cd_rp, a, cd_r1);
    *counter = count_zero_mat_complex_double(cd_r1);
  }

  lu_complex_double(cd_fra, cd_l, cd_u);

  l_step_complex_double(cd_l, cd_frb, cd_y);
  u_step_complex_double(cd_u, cd_y, x);

  for(i=0; i<MATRIX_SIZE; i++){
    xi[i] = x[i];
    xia[i] = x[i];
  }

  iteration_complex_double(cd_fra, cd_l, cd_u, cd_frb, xi);
  iteration_complex_double_another(cd_fr, a, cd_l, cd_u, b, xia);
}

void solve_with_rdft_givens_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE], int *counter){
  int i;
  givens_matrix_list *gml;

  dft_matrix_complex_double(cd_f);
  r_matrix_complex_double(cd_r);
  gml = r_givens_matrix_double();
  mat_mul_complex_double(cd_f,cd_r,cd_fr);
  mat_mul_givens_right_complex_double(cd_fr,gml,cd_frg);
  mat_mul_complex_double(cd_frg,a,cd_fra);
  mat_vec_dot_complex_double(cd_frg,b,cd_frb);

  if(counter != NULL){
    int j;
    for(i=0; i<MATRIX_SIZE; i++){
      for(j=0; j<MATRIX_SIZE; j++)
        cd_r1[i][j] = std::complex<double>(0.0,0.0);
      cd_r1[i][i] = std::complex<double>(1.0,0.0);
    }
    mat_mul_givens_right_complex_double(cd_r1,gml,cd_r2);
    mat_mul_complex_double(cd_r2,a,cd_r3);

    *counter = count_zero_mat_complex_double(cd_r3);
  }

  lu_complex_double(cd_fra, cd_l, cd_u);

  l_step_complex_double(cd_l, cd_frb, cd_y);
  u_step_complex_double(cd_u, cd_y, x);

  for(i=0; i<MATRIX_SIZE; i++){
    xi[i] = x[i];
    xia[i] = x[i];
  }

  iteration_complex_double(cd_fra, cd_l, cd_u, cd_frb, xi);
  iteration_complex_double_another(cd_frg, a, cd_l, cd_u, b, xia);
}
void solve_with_rdft_givens_two_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE], int *counter){
  int i;
  givens_matrix_list *gml1, *gml2;

  dft_matrix_complex_double(cd_f);
  r_matrix_complex_double(cd_r);
  gml1 = r_givens_matrix_double();
  gml2 = r_givens_matrix_double();
  mat_mul_complex_double(cd_f,cd_r,cd_fr);
  mat_mul_givens_right_complex_double(cd_fr,gml1,cd_frg_tmp);
  mat_mul_givens_right_complex_double(cd_frg_tmp,gml2,cd_frg);
  mat_mul_complex_double(cd_frg,a,cd_fra);
  mat_vec_dot_complex_double(cd_frg,b,cd_frb);

  if(counter != NULL){
    int j;
    for(i=0; i<MATRIX_SIZE; i++){
      for(j=0; j<MATRIX_SIZE; j++)
        cd_r1[i][j] = std::complex<double>(0.0,0.0);
      cd_r1[i][i] = std::complex<double>(1.0,0.0);
    }
    mat_mul_givens_right_complex_double(cd_r1,gml1,cd_r2);
    mat_mul_givens_right_complex_double(cd_r2,gml2,cd_r3);
    mat_mul_complex_double(cd_r3,a,cd_r1);
    *counter = count_zero_mat_complex_double(cd_r1);
  }

  lu_complex_double(cd_fra, cd_l, cd_u);

  l_step_complex_double(cd_l, cd_frb, cd_y);
  u_step_complex_double(cd_u, cd_y, x);

  for(i=0; i<MATRIX_SIZE; i++){
    xi[i] = x[i];
    xia[i] = x[i];
  }

  iteration_complex_double(cd_fra, cd_l, cd_u, cd_frb, xi);
  iteration_complex_double_another(cd_frg, a, cd_l, cd_u, b, xia);
}
void iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], givens_matrix_list *gml){
  int itera = 0;
  int iteramax = ITER_MAX;

  while(itera++ < iteramax){
    mat_vec_dot_complex_double(a,x,cd_ax);
    vec_sub_complex_double(b,cd_ax,cd_w);
    l_step_complex_double(l,cd_w,cd_v);
    u_step_complex_double(u,cd_v,cd_z);
    mat_vec_dot_givens_complex_double(gml, cd_z, cd_v);
    vec_add_complex_double(x,cd_v,x);
  }
}
void iteration_both_complex_double_another(std::complex<double> f[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> l[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> u[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], givens_matrix_list *gml){
  int itera = 0;
  int iteramax = ITER_MAX;

  while(itera++ < iteramax){
    mat_vec_dot_complex_double(a,x,cd_ax);
    vec_sub_complex_double(b,cd_ax,cd_w);
    mat_vec_dot_complex_double(f,cd_w,cd_ax);
    l_step_complex_double(l,cd_ax,cd_v);
    u_step_complex_double(u,cd_v,cd_z);
    mat_vec_dot_givens_complex_double(gml, cd_z, cd_v);
    vec_add_complex_double(x,cd_v,x);
  }
}
void solve_with_rdft_both_givens_iteration_complex_double(std::complex<double> a[MATRIX_SIZE][MATRIX_SIZE], std::complex<double> b[MATRIX_SIZE], std::complex<double> x[MATRIX_SIZE], std::complex<double> xi[MATRIX_SIZE], std::complex<double> xia[MATRIX_SIZE], int *counter){
  int i;
  givens_matrix_list *gml1, *gml2;

  dft_matrix_complex_double(cd_f);
  r_matrix_complex_double(cd_r);
  gml1 = r_givens_matrix_double();
  gml2 = r_givens_matrix_double();
  mat_mul_complex_double(cd_f,cd_r,cd_fr);
  mat_mul_givens_right_complex_double(cd_fr, gml1, cd_frg);
  mat_mul_complex_double(cd_frg,a,cd_fra);
  mat_mul_givens_right_complex_double(cd_fra, gml2, cd_tmp);
  mat_vec_dot_complex_double(cd_frg,b,cd_frb);

  if(counter != NULL){
    int j;
    for(i=0; i<MATRIX_SIZE; i++){
      for(j=0; j<MATRIX_SIZE; j++)
        cd_r1[i][j] = std::complex<double>(0.0,0.0);
      cd_r1[i][i] = std::complex<double>(1.0,0.0);
    }
    mat_mul_givens_right_complex_double(cd_r1,gml1,cd_r2);
    mat_mul_complex_double(cd_r2,a,cd_r3);
    mat_mul_givens_right_complex_double(cd_r3,gml2,cd_r1);

    *counter = count_zero_mat_complex_double(cd_r1);
  }

  lu_complex_double(cd_tmp, cd_l, cd_u);

  l_step_complex_double(cd_l, cd_frb, cd_y);
  u_step_complex_double(cd_u, cd_y, x);

  mat_vec_dot_givens_complex_double(gml2,x,cd_y);

  for(i=0; i<MATRIX_SIZE; i++){
    x[i] = cd_y[i];
    xi[i] = x[i];
    xia[i] = x[i];
  }

  iteration_complex_double(cd_fra, cd_l, cd_u, cd_frb, xi, gml2);
  iteration_both_complex_double_another(cd_frg, a, cd_l, cd_u, b, xia, gml2);

  delete_givens_matrix_list(gml1);
  delete_givens_matrix_list(gml2);
}
*/



