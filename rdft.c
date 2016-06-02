#include "rdft.h"
#include "gen.h"
#include <math.h>

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
  char non = 'N', l='L', u='U';

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

  zgemm_(&non,&non, &size,&size,&size, &alpha, atmp,&size, fr,&size, &zero, fra, &size);

  // lu steps
  zgetrfw_(&size, &size, fra, &size);

  // back-forward
  for(int i=0; i<MATRIX_SIZE; i++)
    y[i] = CNUM(b[i],0.0);

  ztrsm_(&l,&l,&non, &u,   &size,&inc, &alpha, fra, &size, y, &size);
  ztrsm_(&l,&u,&non, &non, &size,&inc, &alpha, fra, &size, y, &size);

  // x = Gx
  zgemv_(&non, &size, &size, &alpha, fr, &size, y, &inc, &zero, z, &inc);

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

//-----------------------------------------------
/*
int for_perm[MATRIX_SIZE];
int init_for_perm_flg = 0;

void init_for_perm(){
  int i;
  for(i=0; i<MATRIX_SIZE; i++){
    for_perm[i] = i;
  }
  init_for_perm_flg = 1;
}

void perm_matrix_double(double *r){
  int i,j;
  std::random_device rd;
	std::mt19937 mt(rd());

  if(! init_for_perm_flg){
    init_for_perm();
  }

  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      r[i][j] = 0;
    }
  }

  std::vector<int> perm(for_perm,for_perm+MATRIX_SIZE);
  for(i=0; i<MATRIX_SIZE; i++){
    int idx = mt() % (MATRIX_SIZE - i);
    r[i][perm[idx]] = 1.0;
    perm.erase(perm.begin()+idx);
  }
}
*/
