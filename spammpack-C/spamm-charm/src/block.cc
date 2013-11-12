/** @file
 *
 * The implementation of the SpAMM_Node class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "block.h"
#include "index.h"

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>

/** Update the Frobenius norm. */
void Block::updateNorm (void)
{
  norm = 0;
  norm_2 = 0;
  if(block != NULL)
  {
    for(int i = 0; i < blocksize*blocksize; i++)
    {
      norm_2 += block[i]*block[i];
    }
    norm = sqrt(norm_2);
  }
}

/** The constructor.
 */
Block::Block (void)
{
  blocksize = 0;
  norm = 0;
  norm_2 = 0;
  block = NULL;
}

/** The assignment operator.
 *
 * @param rhs The right hand side of the assignment.
 *
 * @return A new object.
 */
Block & Block::operator= (const Block &rhs)
{
  blocksize = rhs.blocksize;
  block = new double[blocksize*blocksize];
  norm = rhs.norm;
  norm_2 = rhs.norm_2;
  memcpy(block, rhs.block, sizeof(double)*blocksize*blocksize);
}

/** The pup() method.
 *
 * @param p The PUP::er.
 */
void Block::pup (PUP::er &p)
{
  p|blocksize;
  p|norm;
  p|norm_2;
  if(blocksize > 0)
  {
    if(p.isUnpacking())
    {
      block = new double[blocksize*blocksize];
    }
    PUParray(p, block, blocksize*blocksize);
  }
}

/** Set a Block object using a dense matrix.
 *
 * @param blocksize The blocksize.
 * @param A The dense matrix.
 */
void Block::set (const int blocksize, const double *const A)
{
  this->blocksize = blocksize;
  delete[] block;
  block = new double[blocksize*blocksize];
  memcpy(block, A, sizeof(double)*blocksize*blocksize);
  updateNorm();
}

/** Return the blocksize.
 *
 * @return The blocksize.
 */
int Block::getBlocksize (void)
{
  return blocksize;
}

/** Return the square of the norm of this Block.
 *
 * @return The square of the norm.
 */
double Block::getNorm (void)
{
  return norm_2;
}

/** Convert the Block to a dense matrix.
 *
 * @return The dense matrix.
 */
double * Block::toDense (void)
{
  if(blocksize > 0)
  {
    double * dense = new double [blocksize*blocksize];
    memcpy(dense, block, sizeof(double)*blocksize*blocksize);
    return dense;
  }

  else
  {
    return NULL;
  }
}

/** Scale a Block by a factor.
 *
 * @f[ A \leftarrow \alpha A @f]
 *
 * @param alpha The scaling factor.
 */
void Block::scale (const double alpha)
{
  if(blocksize > 0)
  {
    for(int i = 0; i < blocksize*blocksize; i++)
    {
      block[i] *= alpha;
    }
    norm_2 *= alpha*alpha;
    norm = sqrt(norm_2);
  }
}

/** Multiply two Blocks.
 *
 * @f[ A \leftarrow A \times B @f]
 *
 * @param A The matrix A.
 * @param B The matrix B.
 */
void Block::multiply (Block A, Block B)
{
  delete[] block;
  blocksize = A.blocksize;
  block = new double[blocksize*blocksize];
  memset(block, 0, sizeof(double)*blocksize*blocksize);

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

  updateNorm();
}

/** Add a scaled Block.
 *
 * @param alpha The scaling factor.
 * @param A The Block to add.
 */
void Block::add (const double alpha, const double beta, const Block A)
{
  if(blocksize == 0)
  {
    if(A.block != NULL)
    {
      blocksize = A.blocksize;
      block = new double[blocksize*blocksize];
      memset(block, 0, sizeof(double)*blocksize*blocksize);
    }
  }

  if(A.block != NULL)
  {
    for(int i = 0; i < blocksize*blocksize; i++)
    {
      block[i] = alpha*block[i]+beta*A.block[i];
    }
  }

  else
  {
    for(int i = 0; i < blocksize*blocksize; i++)
    {
      block[i] *= alpha;
    }
  }
}

/** Calculate the trace of this Block.
 *
 * @return The trace.
 */
double Block::trace (void)
{
  double trace = 0;

  if(block != NULL)
  {
    for(int i = 0; i < blocksize; i++)
    {
      trace += block[BLOCK_INDEX(i, i, 0, 0, blocksize)];
    }
  }
  return trace;
}

/** Add the identiy matrix to this Block.
 */
void Block::addIdentity (const int blocksize, const double alpha)
{
  if(this->blocksize == 0)
  {
    this->blocksize = blocksize;
    block = new double[blocksize*blocksize];
    memset(block, 0, sizeof(double)*blocksize*blocksize);
  }

  else
  {
    assert(this->blocksize == blocksize);
  }

  for(int i = 0; i < blocksize; i++)
  {
    block[BLOCK_INDEX(i, i, 0, 0, blocksize)] += alpha;
  }

  updateNorm();
}

/** Print a block.
 *
 * @param format The format. Look into printf() as to what to put there.
 */
void Block::print (const char *const format, ...)
{
  va_list ap;
  char buffer[2000];

  va_start(ap, format);
  vsnprintf(buffer, 2000, format, ap);

  printf("%s = [\n", buffer);
  for(int i = 0; i < blocksize; i++) {
    for(int j = 0; j < blocksize; j++)
    {
      printf(" % e", block[BLOCK_INDEX(i, j, 0, 0, blocksize)]);
    }
    printf("\n");
  }
  printf("]\n");
}
