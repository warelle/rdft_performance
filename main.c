#include <stdio.h>
#include <math.h>
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

double *d_x_rdft;
double *d_x_rdft_iter;
double *d_x_rdft_iter_another;
double *d_x_rdft_perm;
double *d_x_rdft_perm_iter;
double *d_x_rdft_perm_iter_another;
double *d_x_rdft_givens;
double *d_x_rdft_givens_iter;
double *d_x_rdft_givens_iter_another;
double *d_x_rdft_givens_two;
double *d_x_rdft_givens_two_iter;
double *d_x_rdft_givens_two_iter_another;
double *d_x_rdft_both_givens;
double *d_x_rdft_both_givens_iter;
double *d_x_rdft_both_givens_iter_another;
double *d_x_gauss;
double *d_x_gauss_iter;
double *d_x_gauss_iter_another;
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
double d_rdft_perm_err;
double d_rdft_perm_iter_err;
double d_rdft_perm_iter_another_err;
double d_rdft_givens_err;
double d_rdft_givens_iter_err;
double d_rdft_givens_iter_another_err;
double d_rdft_givens_two_err;
double d_rdft_givens_two_iter_err;
double d_rdft_givens_two_iter_another_err;
double d_rdft_both_givens_err;
double d_rdft_both_givens_iter_err;
double d_rdft_both_givens_iter_another_err;
double d_gauss_err;
double d_gauss_iter_err;
double d_gauss_iter_another_err;
double d_np_err;
double d_np_iter_err;
double d_np_iter_another_err;
double d_pp_err;
double d_pp_iter_err;
double d_pp_iter_another_err;

void print_double(double d){
  printf("%.10e", d);
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

    alloc_vector_double(&d_x_np, MATRIX_SIZE);

    solve_no_pivoting(d_a, d_b, d_x_np, NULL, NULL);

    daxpy_(&dim, &minus1, d_x, &inc, d_x_np, &inc);
    d_np_err = dnrm2_(&dim, d_x_np, &inc);

    free_vector_double(&d_x_np);
  }
  if(exe & (RDFT | RDFT_ITERATION)){
    int dim = MATRIX_SIZE;
    int inc = 1;
    double minus1 = -1;
    dcomplex *fra=NULL, *frb=NULL, *r=NULL;

    copy_linear_system(d_a,d_x,d_b, a,x,b);

    alloc_vector_double(&d_x_rdft, MATRIX_SIZE);
    alloc_vector_double(&d_x_rdft_iter, MATRIX_SIZE);
    alloc_vector_double(&d_x_rdft_iter_another, MATRIX_SIZE);

    alloc_matrix_complex_double(&fra, MATRIX_SIZE);
    alloc_vector_complex_double(&r, MATRIX_SIZE);
    alloc_vector_complex_double(&frb, MATRIX_SIZE);
    //rdft_original_slow(d_a, d_b, d_x_rdft, d_x_rdft_iter, d_x_rdft_iter_another);
    fftw_rdft_original(d_a, d_b, d_x_rdft, d_x_rdft_iter, d_x_rdft_iter_another, fra, r, frb);

    daxpy_(&dim, &minus1, d_x, &inc, d_x_rdft, &inc);
    daxpy_(&dim, &minus1, d_x, &inc, d_x_rdft_iter, &inc);
    daxpy_(&dim, &minus1, d_x, &inc, d_x_rdft_iter_another, &inc);
    d_rdft_err              = dnrm2_(&dim, d_x_rdft, &inc);
    d_rdft_iter_err         = dnrm2_(&dim, d_x_rdft_iter, &inc);
    d_rdft_iter_another_err = dnrm2_(&dim, d_x_rdft_iter_another, &inc);

    free_matrix_complex_double(&fra);
    free_vector_complex_double(&r);
    free_vector_complex_double(&frb);

    free_vector_double(&d_x_rdft);
    free_vector_double(&d_x_rdft_iter);
    free_vector_double(&d_x_rdft_iter_another);
  }
  /*
  if(exe & (RDFT_PERM | RDFT_PERM_ITERATION)){}
  if(exe & (RDFT_GIVENS | RDFT_GIVENS_ITERATION)){}
  */
  if(exe & (RDFT_GIVENS_TWO | RDFT_GIVENS_TWO_ITERATION)){
    int dim = MATRIX_SIZE;
    int inc = 1;
    double minus1 = -1;
    dcomplex *fra,*r,*frb;
    double *ass;

    copy_linear_system(d_a,d_x,d_b, a,x,b);

    alloc_vector_double(&d_x_rdft_givens_two, MATRIX_SIZE);
    alloc_vector_double(&d_x_rdft_givens_two_iter, MATRIX_SIZE);
    alloc_vector_double(&d_x_rdft_givens_two_iter_another, MATRIX_SIZE);

    alloc_matrix_complex_double(&fra, MATRIX_SIZE);
    alloc_vector_complex_double(&r, MATRIX_SIZE);
    alloc_vector_complex_double(&frb, MATRIX_SIZE);
    alloc_matrix_double(&ass, MATRIX_SIZE);
    fftw_rdft_right_two_givens(d_a, d_b, d_x_rdft_givens_two, d_x_rdft_givens_two_iter, d_x_rdft_givens_two_iter_another, fra,r,frb, ass);

    daxpy_(&dim, &minus1, d_x, &inc, d_x_rdft_givens_two, &inc);
    daxpy_(&dim, &minus1, d_x, &inc, d_x_rdft_givens_two_iter, &inc);
    daxpy_(&dim, &minus1, d_x, &inc, d_x_rdft_givens_two_iter_another, &inc);
    d_rdft_givens_two_err              = dnrm2_(&dim, d_x_rdft_givens_two, &inc);
    d_rdft_givens_two_iter_err         = dnrm2_(&dim, d_x_rdft_givens_two_iter, &inc);
    d_rdft_givens_two_iter_another_err = dnrm2_(&dim, d_x_rdft_givens_two_iter_another, &inc);

    free_matrix_complex_double(&fra);
    free_vector_complex_double(&r);
    free_matrix_complex_double(&frb);
    free_matrix_double(&ass);

    free_vector_double(&d_x_rdft_givens_two);
    free_vector_double(&d_x_rdft_givens_two_iter);
    free_vector_double(&d_x_rdft_givens_two_iter_another);
  }
  /*
  if(exe & (RDFT_BOTH_GIVENS | RDFT_BOTH_GIVENS_ITERATION)){}
  */

  if(exe & (GAUSS | GAUSS_ITERATION)){
    int dim = MATRIX_SIZE;
    int inc = 1;
    double minus1 = -1;
    double *g=NULL,*ga=NULL;

    copy_linear_system(d_a,d_x,d_b, a,x,b);

    alloc_vector_double(&d_x_gauss, MATRIX_SIZE);
    alloc_vector_double(&d_x_gauss_iter, MATRIX_SIZE);
    alloc_vector_double(&d_x_gauss_iter_another, MATRIX_SIZE);

    alloc_matrix_double(&ga, MATRIX_SIZE);
    alloc_matrix_double(&g, MATRIX_SIZE);
    solve_with_gauss_iteration_double(d_a, d_b, d_x_gauss, d_x_gauss_iter, d_x_gauss_iter_another, g,ga);

    daxpy_(&dim, &minus1, d_x, &inc, d_x_gauss, &inc);
    daxpy_(&dim, &minus1, d_x, &inc, d_x_gauss_iter, &inc);
    daxpy_(&dim, &minus1, d_x, &inc, d_x_gauss_iter_another, &inc);
    d_gauss_err              = dnrm2_(&dim, d_x_gauss, &inc);
    d_gauss_iter_err         = dnrm2_(&dim, d_x_gauss_iter, &inc);
    d_gauss_iter_another_err = dnrm2_(&dim, d_x_gauss_iter_another, &inc);

    free_matrix_double(&g);
    free_matrix_double(&ga);

    free_vector_double(&d_x_gauss);
    free_vector_double(&d_x_gauss_iter);
    free_vector_double(&d_x_gauss_iter_another);
  }
  if(exe & (PP | PP_ITERATION)){
    int dim = MATRIX_SIZE;
    int inc = 1;
    double minus1 = -1;

    copy_linear_system(d_a,d_x,d_b, a,x,b);

    alloc_vector_double(&d_x_pp, MATRIX_SIZE);

    solve_with_partial_pivot(d_a, d_b, d_x_pp, NULL, NULL);

    daxpy_(&dim, &minus1, d_x, &inc, d_x_pp, &inc);
    d_pp_err = dnrm2_(&dim, d_x_pp, &inc);

    free_vector_double(&d_x_pp);
  }

  free_matrix_double(&d_a);
  free_vector_double(&d_b);
  free_vector_double(&d_x);
  free_matrix_double(&a);
  free_vector_double(&b);
  free_vector_double(&x);

	// opt
	//   1: graph data
	//   2: readable data
  if(opt == 1){ // graph data
    printf("%d ", x_axis);
    printf("%d ", band_size);
    printf("%d ", 2*band_size-1);
    print_double(d_rdft_err);
    printf(" ");
    print_double(d_rdft_iter_err);
    printf(" ");
    print_double(d_rdft_iter_another_err);
    printf(" ");
    print_double(d_rdft_perm_err);
    printf(" ");
    print_double(d_rdft_perm_iter_err);
    printf(" ");
    print_double(d_rdft_perm_iter_another_err);
    printf(" ");
    print_double(d_rdft_givens_err);
    printf(" ");
    print_double(d_rdft_givens_iter_err);
    printf(" ");
    print_double(d_rdft_givens_iter_another_err);
    printf(" ");
    print_double(d_rdft_givens_two_err);
    printf(" ");
    print_double(d_rdft_givens_two_iter_err);
    printf(" ");
    print_double(d_rdft_givens_two_iter_another_err);
    printf(" ");
    print_double(d_rdft_both_givens_err);
    printf(" ");
    print_double(d_rdft_both_givens_iter_err);
    printf(" ");
    print_double(d_rdft_both_givens_iter_another_err);
    printf(" ");
    print_double(d_gauss_err);
    printf(" ");
    print_double(d_gauss_iter_err);
    printf(" ");
    print_double(d_gauss_iter_another_err);
    printf(" ");
    print_double(d_pp_err);
    printf(" ");
    print_double(d_pp_iter_err);
    printf(" ");
    print_double(d_pp_iter_another_err);
    printf(" ");
    print_double(d_np_err);
    printf("\n");
  }else if(opt == 2){
    if(exe & RDFT){
      printf("RDFT               :");
      print_double(d_rdft_err);
      printf("\n");
    }
    if(exe & RDFT_ITERATION){
      printf("RDFT iteration     :");
      print_double(d_rdft_iter_err);
      printf("\n");
      printf("RDFT iteration     :");
      print_double(d_rdft_iter_another_err);
      printf("\n");
    }
    if(exe & RDFT_PERM){
      printf("RDFT PERM          :");
      print_double(d_rdft_perm_err);
      printf("\n");
    }
    if(exe & RDFT_PERM_ITERATION){
      printf("RDFT PERM iteration:");
      print_double(d_rdft_perm_iter_err);
      printf("\n");
      printf("RDFT PERM iteration:");
      print_double(d_rdft_perm_iter_another_err);
      printf("\n");
    }
    if(exe & RDFT_GIVENS){
      printf("RDFT GVN           :");
      print_double(d_rdft_givens_err);
      printf("\n");
    }
    if(exe & RDFT_GIVENS_ITERATION){
      printf("RDFT GVN iteration :");
      print_double(d_rdft_givens_iter_err);
      printf("\n");
      printf("RDFT GVN iteration :");
      print_double(d_rdft_givens_iter_another_err);
      printf("\n");
    }
    if(exe & RDFT_GIVENS_TWO){
      printf("RDFT GTWO          :");
      print_double(d_rdft_givens_two_err);
      printf("\n");
    }
    if(exe & RDFT_GIVENS_TWO_ITERATION){
      printf("RDFT GTWO iteration:");
      print_double(d_rdft_givens_two_iter_err);
      printf("\n");
      printf("RDFT GTWO iteration:");
      print_double(d_rdft_givens_two_iter_another_err);
      printf("\n");
    }
    if(exe & RDFT_BOTH_GIVENS){
      printf("RDFT BOTH          :");
      print_double(d_rdft_both_givens_err);
      printf("\n");
    }
    if(exe & RDFT_BOTH_GIVENS_ITERATION){
      printf("RDFT BOTH iteration:");
      print_double(d_rdft_both_givens_iter_err);
      printf("\n");
      printf("RDFT BOTH iteration:");
      print_double(d_rdft_both_givens_iter_another_err);
      printf("\n");
    }
    if(exe & GAUSS){
      printf("GAUSS              :");
      print_double(d_gauss_err);
      printf("\n");
    }
    if(exe & GAUSS_ITERATION){
      printf("GAUSS iteration    :");
      print_double(d_gauss_iter_err);
      printf("\n");
      printf("GAUSS iteration    :");
      print_double(d_gauss_iter_another_err);
      printf("\n");
    }
    if(exe & GENP){
      printf("GENP               :");
      print_double(d_np_err);
      printf("\n");
    }
    if(exe & GENP_ITERATION){
      printf("GENP               :");
      print_double(d_np_iter_err);
      printf("\n");
      printf("GENP               :");
      print_double(d_np_iter_another_err);
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
