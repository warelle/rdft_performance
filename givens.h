#ifndef MYGIVENS_H
#define MYGIVENS_H

#include "test.h"
#include "global.h"

//
//   | i  j
// ---------
// i | c -s
// j | s  c
//
typedef struct{
  int i;
  int j;
  double c;
  double s;
} givens_matrix;

typedef struct givens_matrix_list{
  givens_matrix *gm;
  struct givens_matrix_list *next;
}givens_matrix_list;

// THREAD_NUM:2, gm_num:8
//  gsl -> gml -> gm
//   |      |
//   |     gml -> gm
//   |      |
//   |     gml -> gm
//   |      |
//   |     gml -> gm
//   |
//  gsl -> gml -> ...
//         ...
//
typedef struct givens_sequence_list{
  givens_matrix_list *gml;
  struct givens_sequence_list *next;
}givens_sequence_list;


//------------------------------------------------

givens_sequence_list *generate_givens();
void delete_givens(givens_sequence_list *gsl);

//------------------------------------------------

// a = gm * a
void mat_vec_dot_givens_sequence(givens_sequence_list *gsl, double *a);
// A = A * gm
void mat_mul_givens_sequence(givens_sequence_list *gsl, double *a);

#endif
