#ifndef MYGLOBAL_H
#define MYGLOBAL_H

// ----- C99 complex ----- //
#include <complex.h>
typedef double _Complex dcomplex;
#define CNUM(a,b) ((a)+I*(b))

// ----- matrix access ----- //
#define IDX(mat,n, i,j) (mat)[(i)+(j)*(n)]


#endif
