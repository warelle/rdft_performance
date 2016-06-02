*> \brief \b ZGETRF2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       RECURSIVE SUBROUTINE ZGETRFW2( M, N, A, LDA )
*
*       .. Scalar Arguments ..
*       INTEGER            LDA, M, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGETRF2 computes an LU factorization of a general M-by-N matrix A
*> using partial pivoting with row interchanges.
*>
*> The factorization has the form
*>    A = P * L * U
*> where P is a permutation matrix, L is lower triangular with unit
*> diagonal elements (lower trapezoidal if m > n), and U is upper
*> triangular (upper trapezoidal if m < n).
*>
*> This is the recursive version of the algorithm. It divides
*> the matrix into four submatrices:
*>
*>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
*>    A = [ -----|----- ]  with n1 = min(m,n)
*>        [  A21 | A22  ]       n2 = n-n1
*>
*>                                       [ A11 ]
*> The subroutine calls itself to factor [ --- ],
*>                                       [ A12 ]
*>                 [ A12 ]
*> do the swaps on [ --- ], solve A12, update A22,
*>                 [ A22 ]
*>
*> then calls itself to factor A22 and do the swaps on A21.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the M-by-N matrix to be factored.
*>          On exit, the factors L and U from the factorization
*>          A = P*L*U; the unit diagonal elements of L are not stored.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2015
*
*> \ingroup complex16GEcomputational
*
*  =====================================================================
      RECURSIVE SUBROUTINE zgetrfw2( M, N, A, LDA )
*
*  -- LAPACK computational routine (version 3.6.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( lda, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      parameter( one = ( 1.0d+0, 0.0d+0 ),
     $                     zero = ( 0.0d+0, 0.0d+0 ) )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN
      INTEGER            I, N1, N2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            IZAMAX
      EXTERNAL           dlamch, izamax
*     ..
*     .. External Subroutines ..
      EXTERNAL           zgemm, zscal, ztrsm, zerbla
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          max, min
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( m.EQ.0 .OR. n.EQ.0 )
     $   RETURN

      IF( n.EQ.1 ) THEN
*
*        Use unblocked code for one column case
*
*
*        Compute machine safe minimum
*
         sfmin = dlamch('S')
*
*           Compute elements 2:M of the column
*
            IF( abs(a( 1, 1 )) .GE. sfmin ) THEN
               CALL zscal( m-1, one / a( 1, 1 ), a( 2, 1 ), 1 )
            ELSE
               DO 10 i = 1, m-1
                  a( 1+i, 1 ) = a( 1+i, 1 ) / a( 1, 1 )
   10          CONTINUE
            END IF

      ELSE
*
*        Use recursive code
*
         n1 = min( m, n ) / 2
         n2 = n-n1
*
*               [ A11 ]
*        Factor [ --- ]
*               [ A21 ]
*
         CALL zgetrfw2( m, n1, a, lda )
*
*        Solve A12
*
         CALL ztrsm( 'L', 'L', 'N', 'U', n1, n2, one, a, lda,
     $               a( 1, n1+1 ), lda )
*
*        Update A22
*
         CALL zgemm( 'N', 'N', m-n1, n2, n1, -one, a( n1+1, 1 ), lda,
     $               a( 1, n1+1 ), lda, one, a( n1+1, n1+1 ), lda )
*
*        Factor A22
*
         CALL zgetrfw2( m-n1, n2, a( n1+1, n1+1 ), lda )

         END IF
      RETURN
*
*     End of ZGETRF2
*
      END
