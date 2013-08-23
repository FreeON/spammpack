/** @file
 *
 * The blas interface.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __BLAS_INTERFACE_H
#define __BLAS_INTERFACE_H

#include "config.h"

#ifdef DGEMM
/** The dgemm() interface. We are assuming that the dgemm function is using a
 * Fortran style interface, hence all the pointers in the argument list.
 *
 * @param[in] TRANSA
 *           TRANSA is CHARACTER*1
 *            On entry, TRANSA specifies the form of op( A ) to be used in
 *            the matrix multiplication as follows:
 *
 *               TRANSA = 'N' or 'n',  op( A ) = A.
 *
 *               TRANSA = 'T' or 't',  op( A ) = A**T.
 *
 *               TRANSA = 'C' or 'c',  op( A ) = A**T.
 * @param[in]	TRANSB
 *           TRANSB is CHARACTER*1
 *            On entry, TRANSB specifies the form of op( B ) to be used in
 *            the matrix multiplication as follows:
 *
 *               TRANSB = 'N' or 'n',  op( B ) = B.
 *
 *               TRANSB = 'T' or 't',  op( B ) = B**T.
 *
 *               TRANSB = 'C' or 'c',  op( B ) = B**T.
 * @param[in]	M
 *           M is INTEGER
 *            On entry,  M  specifies  the number  of rows  of the  matrix
 *            op( A )  and of the  matrix  C.  M  must  be at least  zero.
 * @param[in]	N
 *           N is INTEGER
 *            On entry,  N  specifies the number  of columns of the matrix
 *            op( B ) and the number of columns of the matrix C. N must be
 *            at least zero.
 * @param[in]	K
 *           K is INTEGER
 *            On entry,  K  specifies  the number of columns of the matrix
 *            op( A ) and the number of rows of the matrix op( B ). K must
 *            be at least  zero.
 * @param[in]	ALPHA
 *           ALPHA is DOUBLE PRECISION.
 *            On entry, ALPHA specifies the scalar alpha.
 * @param[in]	A
 *           A is DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
 *            k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
 *            Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
 *            part of the array  A  must contain the matrix  A,  otherwise
 *            the leading  k by m  part of the array  A  must contain  the
 *            matrix A.
 * @param[in]	LDA
 *           LDA is INTEGER
 *            On entry, LDA specifies the first dimension of A as declared
 *            in the calling (sub) program. When  TRANSA = 'N' or 'n' then
 *            LDA must be at least  max( 1, m ), otherwise  LDA must be at
 *            least  max( 1, k ).
 * @param[in]	B
 *           B is DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
 *            n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
 *            Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
 *            part of the array  B  must contain the matrix  B,  otherwise
 *            the leading  n by k  part of the array  B  must contain  the
 *            matrix B.
 * @param[in]	LDB
 *           LDB is INTEGER
 *            On entry, LDB specifies the first dimension of B as declared
 *            in the calling (sub) program. When  TRANSB = 'N' or 'n' then
 *            LDB must be at least  max( 1, k ), otherwise  LDB must be at
 *            least  max( 1, n ).
 * @param[in]	BETA
 *           BETA is DOUBLE PRECISION.
 *            On entry,  BETA  specifies the scalar  beta.  When  BETA  is
 *            supplied as zero then C need not be set on input.
 * @param[in,out]	C
 *           C is DOUBLE PRECISION array of DIMENSION ( LDC, n ).
 *            Before entry, the leading  m by n  part of the array  C must
 *            contain the matrix  C,  except when  beta  is zero, in which
 *            case C need not be set on entry.
 *            On exit, the array  C  is overwritten by the  m by n  matrix
 *            ( alpha*op( A )*op( B ) + beta*C ).
 * @param[in]	LDC
 *           LDC is INTEGER
 *            On entry, LDC specifies the first dimension of C as declared
 *            in  the  calling  (sub)  program.   LDC  must  be  at  least
 *            max( 1, m ).
 */
extern "C"
void DGEMM (const char *const TRANSA, const char *const TRANSB,
    const int *const M, const int *const N, const int *const K,
    const double *const ALPHA, const double *const A,
    const int *const LDA, const double *const B,
    const int *const LDB, const double *const BETA,
    double *const C, const int *const LDC);
#endif

#endif
