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
extern "C"
void DGEMM (const char *const TRANSA, const char *const TRANSB,
    const int *const M, const int *const N, const int *const K,
    const double *const alpha, const double *const A,
    const int *const LDA, const double *const B,
    const int *const LDB, const double *const beta,
    double *const C, const int *const LDC);
#endif

#endif
