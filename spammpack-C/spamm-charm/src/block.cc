/** @file
 *
 * The implementation of the SpAMM_Node class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "block.h"
#include "index.h"

Block::Block (void)
{
  blocksize = 0;
  block = NULL;
}

void Block::pup (PUP::er &p)
{
  p|blocksize;
}

void Block::set (const int blocksize, const double *const A)
{
  this->blocksize = blocksize;
  delete[] block;
  block = new double[blocksize*blocksize];
  memcpy(block, A, sizeof(double)*blocksize*blocksize);
}

void Block::multiply (Block A, Block B)
{
  delete[] block;
  block = new double[blocksize*blocksize];

#ifdef DGEMM
  double alpha = 1;
  double beta = 1;
  DGEMM("N", "N", &blocksize, &blocksize, &blocksize, &alpha, A.block,
      &blocksize, B.block, &blocksize, &beta, block, &blocksize);
#else
  for(int i = 0; i < blocksize; i++) {
    for(int j = 0; j < blocksize; j++) {
      for(int k = 0; k < blocksize; k++)
      {
        block[BLOCK_INDEX(i, j, 0, 0, blocksize)] +=
          A.block[BLOCK_INDEX(i, k, 0, 0, blocksize)]
          *B.block[BLOCK_INDEX(k, j, 0, 0, blocksize)];
      }
    }
  }
#endif
}
