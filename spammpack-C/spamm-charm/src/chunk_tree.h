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

void *
chunk_tree_alloc (const int N_chunk,
    const int N_basic,
    const int N,
    const int i_lower,
    const int j_lower);

__END_DECLS

#endif
