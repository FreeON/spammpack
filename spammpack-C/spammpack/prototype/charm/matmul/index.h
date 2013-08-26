/** @file
 *
 * Some macros for indexing.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __INDEX_H
#define __INDEX_H

/** The linear offset inside a dense square matrix block. Column-major order. */
#define BLOCK_INDEX(i, j, iLower, jLower, blocksize) ((i-iLower)+(j-jLower)*blocksize)

/** The linear tree index. */
#define CHILD_INDEX(i, j) ((i << 1) | j)

#endif
