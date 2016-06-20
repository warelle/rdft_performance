#include "matrix_complex.h"
#include "rdft.h"
#include "givens.h"
#include "gen.h"
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <stdio.h>

#define MATH_PI 3.14159265358979

//-----------------------------------------------
void dft_matrix_complex_double(dcomplex *f){
  int i, j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      IDX(f,MATRIX_SIZE,i,j) = CNUM(cos(-2.0*MATH_PI*i*j/MATRIX_SIZE), sin(-2.0*MATH_PI*i*j/MATRIX_SIZE));
    }
  }
}

void r_matrix_complex_double(dcomplex *r){
  int i;
  for(i=0; i<MATRIX_SIZE; i++){
    double rval = uniform()*MATH_PI*(rand_normal_double(0,1) > 0 ? -1 : 1);
    r[i] = CNUM(cos(-2.0*MATH_PI*rval/MATRIX_SIZE), sin(-2.0*MATH_PI*rval/MATRIX_SIZE) );
  }
}
//-----------------------------------------------
void rdft_original_slow(double *a, double *b, double *x, double *xi, double *xia){
  int size=MATRIX_SIZE, inc=1;
  dcomplex *fra, *fr, *r, *y, *z, *atmp;
  dcomplex alpha=CNUM(1.0, 0.0), zero=CNUM(0.0, 0.0);

  alloc_matrix_complex_double(&fra, MATRIX_SIZE);
  alloc_matrix_complex_double(&fr, MATRIX_SIZE);
  alloc_vector_complex_double(&r, MATRIX_SIZE);
  alloc_vector_complex_double(&y, MATRIX_SIZE);
  alloc_vector_complex_double(&z, MATRIX_SIZE);
  alloc_matrix_complex_double(&atmp, MATRIX_SIZE);

  dft_matrix_complex_double(fr);
  r_matrix_complex_double(r);

  for(int i=0; i<MATRIX_SIZE; i++)
    for(int j=0; j<MATRIX_SIZE; j++)
      IDX(atmp, MATRIX_SIZE, i,j) = CNUM(IDX(a, MATRIX_SIZE, i,j),0.0);

  for(int i=0; i<MATRIX_SIZE; i++)
    zscal_(&size, (r+i), fr, &inc);

  zgemm_("N","N", &size,&size,&size, &alpha, atmp,&size, fr,&size, &zero, fra, &size);

  // lu steps
  zgetrfw_(&size, &size, fra, &size);

  // back-forward
  for(int i=0; i<MATRIX_SIZE; i++)
    y[i] = CNUM(b[i],0.0);

  ztrsm_("L","L","N", "U", &size,&inc, &alpha, fra, &size, y, &size);
  ztrsm_("L","U","N", "N", &size,&inc, &alpha, fra, &size, y, &size);

  // x = Gx
  zgemv_("N", &size, &size, &alpha, fr, &size, y, &inc, &zero, z, &inc);

  for(int i=0; i<MATRIX_SIZE; i++)
    x[i] = creal(z[i]);

  // iteration_double(d_fa, d_l, d_u, d_fb, xi);
  // iteration_double_another(d_f, a, d_l, d_u, b, xia);

  free_matrix_complex_double(&fra);
  free_matrix_complex_double(&fr);
  free_vector_complex_double(&r);
  free_vector_complex_double(&y);
  free_vector_complex_double(&z);
  free_matrix_complex_double(&atmp);
}

void fftw_rdft_original(double *a, double *b, double *x, double *xi, double *xia, dcomplex *fra, dcomplex *r, dcomplex *frb){
  /*
   * FRA = RFA =R(FA)
   * FRAx = FRb <=> R(FA)x = R(Fb)
   *
   */
  int size=MATRIX_SIZE, inc=1;
  dcomplex alpha=CNUM(1.0, 0.0);

  // FFT
  //#pragma omp parallel for
  for(int i=0; i<MATRIX_SIZE; i++){
    fftw_plan ftplan = fftw_plan_dft_r2c_1d(MATRIX_SIZE, a+(MATRIX_SIZE*i), fra+(MATRIX_SIZE*i), FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    fftw_execute(ftplan);
    fftw_destroy_plan(ftplan);

    for(int j=MATRIX_SIZE-1; j>=MATRIX_SIZE/2; j--)
      fra[MATRIX_SIZE*i+j] = conj(fra[MATRIX_SIZE*i + MATRIX_SIZE-j]);
  }
  fftw_plan ftplan_b = fftw_plan_dft_r2c_1d(MATRIX_SIZE, b, frb, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  fftw_execute(ftplan_b);
  fftw_destroy_plan(ftplan_b);
  for(int j=MATRIX_SIZE-1; j>=MATRIX_SIZE/2; j--)
    frb[j] = conj(frb[MATRIX_SIZE-j]);

  // Multiply R=diag(...)
  r_matrix_complex_double(r);
  #pragma omp parallel for
  for(int i=0; i<MATRIX_SIZE; i++)
    for(int j=0; j<MATRIX_SIZE; j++)
      fra[i*MATRIX_SIZE + j] *= r[j];
  for(int i=0; i<MATRIX_SIZE; i++)
    frb[i] = r[i]*frb[i];

  // lu steps
  zgetrfw_(&size, &size, fra, &size);

  ztrsm_("L","L","N", "U", &size,&inc, &alpha, fra, &size, frb, &size);
  ztrsm_("L","U","N", "N", &size,&inc, &alpha, fra, &size, frb, &size);

  for(int i=0; i<MATRIX_SIZE; i++)
    x[i] = creal(frb[i]);

  // iteration_double(d_fa, d_l, d_u, d_fb, xi);
  // iteration_double_another(d_f, a, d_l, d_u, b, xia);
}

void fftw_rdft_right_two_givens(double *a, double *b, double *x, double *xi, double *xia, dcomplex *fra, dcomplex *r, dcomplex *frb, double *ass){
  /*
   * FRASS = RFASS =R(F(ASS))
   * FRASSy = FRb <=> R(F(ASS)y) = R(Fb), x=SSy
   *
   */
  int size=MATRIX_SIZE, size2=MATRIX_SIZE*MATRIX_SIZE, inc=1;
  dcomplex alpha=CNUM(1.0, 0.0);

  givens_sequence_list *gsl1 = generate_givens(),
                       *gsl2 = generate_givens();

  dcopy_(&size2, a, &inc, ass, &inc);
  mat_mul_givens_sequence(gsl1, ass);
  mat_mul_givens_sequence(gsl2, ass);

  // FFT
  //#pragma omp parallel for
  for(int i=0; i<MATRIX_SIZE; i++){
    fftw_plan ftplan = fftw_plan_dft_r2c_1d(MATRIX_SIZE, ass+(MATRIX_SIZE*i), fra+(MATRIX_SIZE*i), FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    fftw_execute(ftplan);
    fftw_destroy_plan(ftplan);

    for(int j=MATRIX_SIZE-1; j>=MATRIX_SIZE/2; j--)
      fra[MATRIX_SIZE*i+j] = conj(fra[MATRIX_SIZE*i + MATRIX_SIZE-j]);
  }
  fftw_plan ftplan_b = fftw_plan_dft_r2c_1d(MATRIX_SIZE, b, frb, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  fftw_execute(ftplan_b);
  fftw_destroy_plan(ftplan_b);
  for(int j=MATRIX_SIZE-1; j>=MATRIX_SIZE/2; j--)
    frb[j] = conj(frb[MATRIX_SIZE-j]);

  // Multiply R=diag(...)
  r_matrix_complex_double(r);
  #pragma omp parallel for
  for(int i=0; i<MATRIX_SIZE; i++)
    for(int j=0; j<MATRIX_SIZE; j++)
      fra[i*MATRIX_SIZE + j] *= r[j];
  for(int i=0; i<MATRIX_SIZE; i++)
    frb[i] = r[i]*frb[i];

  // lu steps
  zgetrfw_(&size, &size, fra, &size);

  ztrsm_("L","L","N", "U", &size,&inc, &alpha, fra, &size, frb, &size);
  ztrsm_("L","U","N", "N", &size,&inc, &alpha, fra, &size, frb, &size);

  for(int i=0; i<MATRIX_SIZE; i++)
    x[i] = creal(frb[i]);

  mat_vec_dot_givens_sequence(gsl2,x);
  mat_vec_dot_givens_sequence(gsl1,x);

  // iteration_double(d_fa, d_l, d_u, d_fb, xi);
  // iteration_double_another(d_f, a, d_l, d_u, b, xia);
}

