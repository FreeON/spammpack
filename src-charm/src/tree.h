/** @file
 *
 * The header file for parallel trees.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __TREE_H
#define __TREE_H

#include <stdlib.h>

__BEGIN_DECLS

void *
tree_new (const int N, const int N_chunk);

__END_DECLS

#endif
