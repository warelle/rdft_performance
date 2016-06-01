#ifndef MYGEN_H
#define MYGEN_H

#include "test.h"

// return: random value
double uniform();
double rand_normal_double(double m, double s);

// return: mat,vec
void generate_linear_system(double *a, double *x, double *b, double range, int band_size);

// matrices
void gen_identity_float128(__float128 *mat);
void gen_random_float128(__float128 *mat, __float128 range);
void gen_arrowhead_float128(__float128 *mat, __float128 range);
void gen_diag_big_float128(__float128 *mat, __float128 range);
void gen_band_no_side_float128(__float128 *mat, __float128 range, int band_size);
void gen_band_float128(__float128 *mat, __float128 range, int band_size);
void gen_normal_float128(__float128 *mat);
void gen_uniform_absone_float128(__float128 *mat);
void gen_uniform_positive_float128(__float128 *mat);
void gen_set_absone_float128(__float128 *mat);
void gen_set_abspositive_float128(__float128 *mat);
void gen_ijabs_float128(__float128 *mat);
void gen_maxij_float128(__float128 *mat);
void gen_gfpp_float128(__float128 *mat);
void gen_fiedler_float128(__float128 *mat);
void gen_hadamard_float128(__float128 *mat);

#endif
