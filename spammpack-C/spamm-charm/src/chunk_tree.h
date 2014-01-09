/** @file
 *
 * The header file for chunk functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __CHUNK_TREE_H
#define __CHUNK_TREE_H

#include <stdlib.h>

__BEGIN_DECLS

size_t
chunk_tree_sizeof (const int N_chunk, const int N_basic);

void *
chunk_tree_alloc (const int N_chunk,
    const int N_basic,
    const int N,
    const int i_lower,
    const int j_lower);

void
chunk_tree_set (void *const chunk, const double *const A);

double
chunk_tree_get_norm (const void *const chunk);

void
chunk_tree_multiply (const double tolerance,
    const void *const A,
    const void *const B,
    void *const C);

double *
chunk_tree_to_dense (const void *const chunk);

__END_DECLS

#endif
