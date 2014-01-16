/** @file
 *
 * The header file for chunk functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __CHUNK_BLOCK_H
#define __CHUNK_BLOCK_H

#include <stdlib.h>

__BEGIN_DECLS

void
chunk_block_multiply (const double *const restrict A,
    const double *const restrict B,
    double *const restrict C,
    const int N);

__END_DECLS

#endif
