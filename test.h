#ifndef MYTEST_H
#define MYTEST_H

#define MATRIX_SIZE 4096
#define ITER_MAX 1
//#define THDEAD_NUM 2

#include "global.h"
#include "blas_extern.h"

// -------------------------------------- //
#define NONE                        0
#define RDFT                        1
#define RDFT_ITERATION              2
#define GENP                        4
#define GENP_ITERATION              8
#define GAUSS                      16
#define GAUSS_ITERATION            32
#define PP                         64
#define PP_ITERATION               128
#define RDFT_PERM                  256
#define RDFT_PERM_ITERATION        512
#define RDFT_GIVENS                1024
#define RDFT_GIVENS_ITERATION      2048
#define RDFT_GIVENS_TWO            4096
#define RDFT_GIVENS_TWO_ITERATION  8192
#define RDFT_BOTH_GIVENS           16384
#define RDFT_BOTH_GIVENS_ITERATION 32768
#define ALL                        (RDFT | RDFT_ITERATION | RDFT_PERM | RDFT_PERM_ITERATION | RDFT_GIVENS | RDFT_GIVENS_ITERATION | GAUSS | GAUSS_ITERATION | PP | PP_ITERATION | GENP | GENP_ITERATION)

// [1] EXECUTION METHOD
#define EXE (GENP | PP | GAUSS | RDFT | RDFT_GIVENS_TWO)
// -------------------------------------- //

#define RANDOM        1
#define DIAG_BIG     (1 << 1)
#define ARROWHEAD    (1 << 2)
#define IDENTITY     (1 << 3)
#define BAND_NO_SIDE (1 << 4)
#define BAND         (1 << 5)

#define NORMAL           (1 << 6)
#define UNIFORM_ABSONE   (1 << 7)
#define UNIFORM_POSITIVE (1 << 8)
#define SET_ABSONE       (1 << 9)
#define SET_ABSPOSITIVE  (1 << 10)

#define IJABS (1 << 11)
#define MAXIJ (1 << 12)

#define GFPP     (1 << 13)
#define FIEDLER  (1 << 14)
#define HADAMARD (1 << 15)

// [2] TEST MATRIX
#define TEST_MATRIX RANDOM
// -------------------------------------- //

#ifdef THREAD_NUM
#include <omp.h>
omp_set_num_threads(THREAD_NUM);
#endif

#endif
