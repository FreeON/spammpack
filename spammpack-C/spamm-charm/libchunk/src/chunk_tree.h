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

int
chunk_tree_get_N_chunk (void *const chunk);

int
chunk_tree_get_N_basic (void *const chunk);

double
chunk_tree_get_norm_2 (const void *const chunk);

void
chunk_tree_print (const void *const chunk,
    const char *const format, ...);

void
chunk_tree_add (const double alpha, void *const A,
    const double beta, const void *const B);

void
chunk_tree_multiply (const double tolerance,
    const void *const A,
    const void *const B,
    void *const C,
    const short tree_only);

double
chunk_tree_trace (const void *const chunk);

void
chunk_tree_scale (const double alpha, void *const chunk);

void
chunk_tree_add_identity (const double alpha, void *const chunk);

double *
chunk_tree_to_dense (const void *const chunk);

void
chunk_tree_delete (void **const chunk);

size_t
chunk_tree_get_complexity (const void *const chunk);

__END_DECLS

#endif
