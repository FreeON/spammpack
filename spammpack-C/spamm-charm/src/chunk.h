/** @file
 *
 * The header file for chunk functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __CHUNK_H
#define __CHUNK_H

#include <stdlib.h>

#ifdef __cplusplus
#define __BEGIN_DECLARATIONS extern "C" {
#define __END_DECLARATIONS }
#else
#define __BEGIN_DECLARATIONS
#define __END_DECLARATIONS
#endif

__BEGIN_DECLARATIONS

size_t
chunk_sizeof (const int N_chunk, const int N_basic);

void *
chunk_alloc (const int N_chunk,
    const int N_basic,
    const int N,
    const unsigned int i_lower,
    const unsigned int j_lower);

void
chunk_set (void *const chunk, const double *const A);

double
chunk_get_norm (const void *const chunk);

void
chunk_print (const void *const chunk,
    const char *const format, ...);

void
chunk_add (const double alpha, void *const A,
    const double beta, const void *const B);

void
chunk_multiply (const void *const A, const void *const B, void *const C);

double
chunk_trace (const void *const chunk);

void
chunk_scale (const double alpha, void *const chunk);

void
chunk_add_identity (const double alpha, void *const chunk);

double *
chunk_to_dense (const void *const chunk);

__END_DECLARATIONS

#endif
