#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "matrix_double.h"
#include "matrix_complex.h"
#include "test.h"
#include "lu.h"
#include "gen.h"
#include "solve.h"
#include "rdft.h"
#include "gauss.h"

// double
double *a;
double *b;
double *x;

dcomplex *dc_aa;
dcomplex *dc_bb;
dcomplex *dc_xx;
dcomplex *dc_xi;
dcomplex *dc_xa;
double *d_a;
double *d_b;
double *d_x;

double *d_x_pp;
double *d_x_pp_iter;
double *d_x_pp_iter_another;
double *d_x_pp_with_round_error;
double *d_x_pp_iter_with_round_error;
double *d_x_pp_iter_another_with_round_error;
double *d_x_np;
double *d_x_np_iter;
double *d_x_np_iter_another;

double d_rdft_err;
double d_rdft_iter_err;
double d_rdft_iter_another_err;
double d_np_err;
double d_np_iter_err;
double d_np_iter_another_err;
double d_pp_err;
double d_pp_iter_err;
double d_pp_iter_another_err;

void print_double(double d){
  printf("%.10e", d);
}

void solve_np(double *a, double *b){
  double alpha=1.0;
  int size = MATRIX_SIZE, inc=1;
  char non = 'N', l='L', u='U';

  // lu steps
  dgetrfw_(&size, &size, a, &size);

  dtrsm_(&l,&l,&non, &u,   &size,&inc, &alpha, a, &size, b, &size);
  dtrsm_(&l,&u,&non, &non, &size,&inc, &alpha, a, &size, b, &size);
}
void solve_znp(dcomplex *a, dcomplex *b){
  int size=MATRIX_SIZE, inc=1;
  dcomplex alpha=CNUM(1.0, 0.0);
  char non = 'N', l='L', u='U';

  // lu steps
  zgetrfw_(&size, &size, a, &size);

  // back-forward
  ztrsm_(&l,&l,&non, &u,   &size,&inc, &alpha, a, &size, b, &size);
  ztrsm_(&l,&u,&non, &non, &size,&inc, &alpha, a, &size, b, &size);

}
void solve_pp(double *aa, double *bb, double *xx, int *piv, double *work, float *workf){
  int dim = MATRIX_SIZE;
  int nrhs = 1;
  int iter,info;

  dsgesv_(&dim, &nrhs, aa,&dim, piv, bb,&dim, xx,&dim, work, workf, &iter, &info);
}

void run(int dat, int opt, int exe, int band_size, int x_axis){
  alloc_matrix_double(&a,MATRIX_SIZE);
  alloc_vector_double(&b,MATRIX_SIZE);
  alloc_vector_double(&x,MATRIX_SIZE);

  alloc_matrix_double(&d_a,MATRIX_SIZE);
  alloc_vector_double(&d_b,MATRIX_SIZE);
  alloc_vector_double(&d_x,MATRIX_SIZE);

  generate_linear_system(a,x,b, 1.0, band_size);

  if(exe & (GENP | GENP_ITERATION)){
    int dim = MATRIX_SIZE;
    int inc = 1;
    double minus1 = -1;

    copy_linear_system(d_a,d_x,d_b, a,x,b);

    solve_np(d_a, d_b);

    daxpy_(&dim, &minus1, d_x, &inc, d_b, &inc);
    d_np_err = dnrm2_(&dim, d_b, &inc);
  }
  if(exe & (PP | PP_ITERATION)){
    int dim = MATRIX_SIZE;
    int inc = 1;
    double minus1 = -1;

    copy_linear_system(d_a,d_x,d_b, a,x,b);

    int *piv = (int*)malloc(sizeof(int)*MATRIX_SIZE);
    double *work  = (double*)malloc(sizeof(double)*MATRIX_SIZE);
    float  *workf = (float*) malloc(sizeof(float)*MATRIX_SIZE*(MATRIX_SIZE+1));
    alloc_vector_double(&d_x_pp,MATRIX_SIZE);

    solve_pp(d_a, d_b, d_x_pp, piv, work, workf);

    daxpy_(&dim, &minus1, d_x, &inc, d_x_pp, &inc);
    d_pp_err = dnrm2_(&dim, d_x_pp, &inc);
    free_vector_double(&d_x_pp);
  }
  if(exe & (GENP | GENP_ITERATION)){
    int dim = MATRIX_SIZE;
    int inc = 1;
    double minus1 = -1;
    double *xx;
    dcomplex *atmp, *xtmp;
    alloc_matrix_complex_double(&atmp, MATRIX_SIZE);
    alloc_vector_complex_double(&xtmp, MATRIX_SIZE);
    alloc_vector_double(&xx, MATRIX_SIZE);

    copy_linear_system(d_a,d_x,d_b, a,x,b);

    for(int i=0; i<MATRIX_SIZE; i++)
      for(int j=0; j<MATRIX_SIZE; j++)
        IDX(atmp, MATRIX_SIZE, i,j) = CNUM(IDX(a, MATRIX_SIZE, i,j),0.0);
    for(int i=0; i<MATRIX_SIZE; i++)
      xtmp[i] = CNUM(d_b[i],0.0);

    solve_znp(atmp, xtmp);

    for(int i=0; i<MATRIX_SIZE; i++)
      xx[i] = creal(xtmp[i]);

    daxpy_(&dim, &minus1, d_x, &inc, xx, &inc);
    d_rdft_err = dnrm2_(&dim, xx, &inc);

    free_vector_double(&xx);
    free_vector_complex_double(&xtmp);
    free_matrix_complex_double(&atmp);
  }


  free_matrix_double(&d_a);
  free_vector_double(&d_b);
  free_vector_double(&d_x);
  free_matrix_double(&a);
  free_vector_double(&b);
  free_vector_double(&x);

  if(opt == 2){
    if(exe & RDFT){
      printf("RDFT               :");
      print_double(d_rdft_err);
      printf("\n");
    }
    if(exe & GENP){
      printf("GENP               :");
      print_double(d_np_err);
      printf("\n");
    }
    if(exe & PP){
      printf("Partial Pivot      :");
      print_double(d_pp_err);
      printf("\n");
    }
  }
}

int main(){
  //for(int i=0; i<MATRIX_SIZE/2; i++){
  for(int i=0; i<1; i++){
//    for(int j=0; j<3; j++)
      run(i,2, EXE, i+1, MATRIX_SIZE);
    //fprintf(stderr, "%d ", i);
    //fprintf(stderr, "%d ", MATRIX_SIZE);
  }

  return 0;
}
