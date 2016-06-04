#include "gen.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <quadmath.h>

// --------------------------------------------------------
int *ccc;
static int ccc_flag = 0;

void ccc_init(){
  ccc_flag = 1;
  ccc = (int*)malloc(sizeof(int)*(MATRIX_SIZE+2)*(MATRIX_SIZE+2));
}
void ccc_fin(){
  ccc_flag = 0;
  free(ccc);
}
int combination(int i, int j){
  if(i < 0 || j < 0){
    return 1;
  }
  if(IDX(ccc,MATRIX_SIZE,i,j) != 0){
    return IDX(ccc,MATRIX_SIZE,i,j);
  }
  if(j == 0 || i <= j){
    IDX(ccc,MATRIX_SIZE,i,j) = 1;
    return IDX(ccc,MATRIX_SIZE,i,j);
  }
  IDX(ccc,MATRIX_SIZE,i,j) = combination(i-1,j-1) + combination(i-1,j);
  if(i > j)
    IDX(ccc,MATRIX_SIZE,i,i-j) = IDX(ccc,MATRIX_SIZE,i,j);
  return IDX(ccc,MATRIX_SIZE,i,j);
}
unsigned msb(unsigned n){
  n|=n>>1; n|=n>>2; n|=n>>4;
  n|=n>>8; n|=n>>16;
  return n^(n>>1);
}
int hadamard(int y,int x){
  int m,n;
  if(x==0||y==0)return  1;
  if(x==1&&y==1)return -1;
  n=(m=msb(x|y))-1;
  m=(m<=x&&m<=y)?-1:1;
  return m*hadamard(y&n,x&n);
}
// --------------------------------------------------------

// return: random value
double uniform(){
  static int init_flg = 0;
  if(!init_flg){
    init_flg = 1;
    srand((unsigned)time(NULL));
  }
  return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}
double rand_normal_double(double m, double s){
    double z = sqrt(-2.0*log(uniform()))*sin(2.0*M_PI*uniform());
    return m + s*z;
}
__float128 rand_normal_float_128(__float128 m, __float128 s){
    __float128 z = sqrtq(-2.0*logq(uniform()))*sinq(2.0*M_PIq*uniform());
    return m + s*z;
}
__float128 random_float_128(__float128 range){
  __float128 r = 0.0Q;
  r += uniform();
  r *= 1.0e-10;
  r += uniform();
  r *= 1.0e-10;
  r += uniform();
  r *= range;
  return uniform() > 0.5 ? r : -r;
}

// --------------------------------------------------------
// return: vec
void generate_vector_float128(__float128 *vec, __float128 range){
  int i;
  for(i=0; i<MATRIX_SIZE; i++)
    vec[i] = random_float_128(range);
}

// return: mat
void generate_matrix_float128(__float128 *mat, __float128 range, int band_size){
  gen_random_float128(mat,range);
//  gen_diag_big_float128(mat,range);
//  gen_arrowhead_float128(mat,range);
//  gen_identity_float128(mat);
//  gen_band_no_side_float128(mat, range, opt);
//  gen_band_float128(mat, range, band_size);

//  gen_normal_float128(mat);
//  gen_uniform_absone_float128(mat);
//  gen_uniform_positive_float128(mat);
//  gen_set_absone_float128(mat);
//  gen_set_abspositive_float128(mat);

//  gen_ijabs_float128(mat);
//  gen_maxij_float128(mat);

//  gen_gfpp_float128(mat);
//  gen_fiedler_float128(mat);
//  gen_hadamard_float128(mat);

}

// return: mat,vec
void generate_linear_system_float_128(__float128 *mat, __float128 *vec, __float128 range, int band_size){
  generate_matrix_float128(mat, range, band_size);
  generate_vector_float128(vec, range);
}
void generate_linear_system(double *a, double *x, double *b, double range, int band_size){
  __float128 *mat = (__float128*)malloc((sizeof(__float128))*MATRIX_SIZE*MATRIX_SIZE);
  __float128 *vec = (__float128*)malloc((sizeof(__float128))*MATRIX_SIZE);
  __float128 *mv  = (__float128*)malloc((sizeof(__float128))*MATRIX_SIZE);

  if(mat == NULL || vec == NULL || mv == NULL){
    fprintf(stderr,"malloc error in gen.c");
  }

  generate_linear_system_float_128(mat, vec, (__float128)range, band_size);

  for(int i=0; i<MATRIX_SIZE; i++){
    mv[i] = 0.0;
    for(int j=0; j<MATRIX_SIZE; j++){
      mv[i] = mv[i] + IDX(mat,MATRIX_SIZE,i,j)*vec[j];
    }
  }

  for(int i=0; i<MATRIX_SIZE; i++){
    for(int j=0; j<MATRIX_SIZE; j++){
      IDX(a,MATRIX_SIZE,i,j) = (double)IDX(mat,MATRIX_SIZE,i,j);
    }
    b[i] = (double)mv[i];
    x[i] = (double)vec[i];
  }

  free(mat);
  free(vec);
  free(mv);
}

void copy_linear_system(double *d_a, double *d_x, double *d_b, double *a, double *x, double *b){
  for(int i=0; i<MATRIX_SIZE; i++){
    for(int j=0; j<MATRIX_SIZE; j++){
      IDX(d_a,MATRIX_SIZE, i,j) = IDX(a,MATRIX_SIZE,i,j);
    }
    d_x[i] = x[i];
    d_b[i] = b[i];
  }
}


// --------------------------------------------------------

/*--- specific type of matrices ---*/
void gen_identity_float128(__float128 *mat){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      if(i == j)
        IDX(mat,MATRIX_SIZE,i,j) = 1.0Q;
      else
        IDX(mat,MATRIX_SIZE,i,j) = 0.0Q;
}
void gen_random_float128(__float128 *mat, __float128 range){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      IDX(mat,MATRIX_SIZE,i,j) = random_float_128(range);
}
void gen_arrowhead_float128(__float128 *mat, __float128 range){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      if(i == j || i==0 || j==0)
        IDX(mat,MATRIX_SIZE,i,j) = random_float_128(range);
      else
        IDX(mat,MATRIX_SIZE,i,j) = 0.0Q;
}
void gen_diag_big_float128(__float128 *mat, __float128 range){
  int i,j;
  __float128 big_f = 1000000000000000000.0Q;
  __float128 small_f = 1.0Q;
  for(i=0; i<MATRIX_SIZE; i++)
    for(j=0; j<MATRIX_SIZE; j++)
      if(i == j)
        IDX(mat,MATRIX_SIZE,i,j) = random_float_128(range*big_f);
      else
        IDX(mat,MATRIX_SIZE,i,j) = random_float_128(range*small_f);
}
// [NOTE] BAND_SIZE will be 2*band_size-1
void gen_band_no_side_float128(__float128 *mat, __float128 range, int band_size){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    if(i<band_size){
      IDX(mat,MATRIX_SIZE,i,0) = random_float_128(range);
      IDX(mat,MATRIX_SIZE,0,i) = random_float_128(range);
    }else{
      IDX(mat,MATRIX_SIZE,i,0) = 0.0Q;
      IDX(mat,MATRIX_SIZE,0,i) = 0.0Q;
    }
  }
  for(i=1; i<MATRIX_SIZE; i++)
    for(j=1; j<MATRIX_SIZE; j++)
      if(IDX(mat,MATRIX_SIZE,i-1,j-1) != 0.0Q)
        IDX(mat,MATRIX_SIZE,i,j) = random_float_128(range);
      else
        IDX(mat,MATRIX_SIZE,i,j) = 0.0Q;
}
void gen_band_float128(__float128 *mat, __float128 range, int band_size){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    if(i<band_size || (MATRIX_SIZE-band_size) < i){
      IDX(mat,MATRIX_SIZE,i,0) = random_float_128(range);
      IDX(mat,MATRIX_SIZE,0,i) = random_float_128(range);
    }else{
      IDX(mat,MATRIX_SIZE,i,0) = 0.0Q;
      IDX(mat,MATRIX_SIZE,0,i) = 0.0Q;
    }
  }
  for(i=1; i<MATRIX_SIZE; i++)
    for(j=1; j<MATRIX_SIZE; j++)
      if(IDX(mat,MATRIX_SIZE,i-1,j-1) != 0.0Q)
        IDX(mat,MATRIX_SIZE,i,j) = random_float_128(range);
      else
        IDX(mat,MATRIX_SIZE,i,j) = 0.0Q;
}


/*--- specific type of matrices ---*/
void gen_normal_float128(__float128 *mat){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      IDX(mat,MATRIX_SIZE,i,j) = rand_normal_double(0,1);
    }
  }
}
void gen_uniform_absone_float128(__float128 *mat){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      IDX(mat,MATRIX_SIZE,i,j) = uniform() * (rand_normal_double(0,1) > 0 ? -1 : 1);
    }
  }
}
void gen_uniform_positive_float128(__float128 *mat){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      IDX(mat,MATRIX_SIZE,i,j) = uniform();
    }
  }
}
void gen_set_absone_float128(__float128 *mat){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      IDX(mat,MATRIX_SIZE,i,j) = rand_normal_double(0,1) > 0 ? -1 : 1;
    }
  }
}
void gen_set_abspositive_float128(__float128 *mat){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      IDX(mat,MATRIX_SIZE,i,j) = rand_normal_double(0,1) > 0 ? 0 : 1;
    }
  }
}

void gen_ijabs_float128(__float128 *mat){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      IDX(mat,MATRIX_SIZE,i,j) = (i-j)>0 ? i-j:j-i;
    }
  }
}
void gen_maxij_float128(__float128 *mat){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      IDX(mat,MATRIX_SIZE,i,j) = i>j ? i : j;
    }
  }
}

void gen_gfpp_float128(__float128 *mat){
  int i,j;
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      if(i<j){
        IDX(mat,MATRIX_SIZE,i,j) = 0;
      }else if(i==j){
        IDX(mat,MATRIX_SIZE,i,j) = 1.0;
      }else{
        IDX(mat,MATRIX_SIZE,i,j) = -1.0;
      }
    }
    IDX(mat,MATRIX_SIZE,i,MATRIX_SIZE-1) = 1.0;
  }
}

void gen_fiedler_float128(__float128 *mat){
  int i,j,k;
  for(i=0; i<MATRIX_SIZE; i++){
    int ret_flg = 0;
    k = i;
    for(j=0; j<MATRIX_SIZE; j++){
      IDX(mat,MATRIX_SIZE,i,j) = k;
      if(ret_flg == 0){
        k--;
        if(k < 0){
          ret_flg = 1;
          k = 0;
        }
      }else{
        k++;
      }
    }
  }
}

void gen_hadamard_float128(__float128 *mat){
  int i,j;
  ccc_init();
  for(i=0; i<MATRIX_SIZE; i++){
    for(j=0; j<MATRIX_SIZE; j++){
      IDX(mat,MATRIX_SIZE,i,j) = hadamard(i,j);
    }
  }
  ccc_fin();
}


