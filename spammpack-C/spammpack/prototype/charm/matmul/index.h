#ifndef __INDEX_H
#define __INDEX_H

#define BLOCK_INDEX(i, j, iLower, jLower, blocksize) ((i-iLower)*blocksize+(j-jLower))
#define CHILD_INDEX(i, j) ((i << 1) | j)

#endif
