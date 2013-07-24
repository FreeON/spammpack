/** @file
 *
 * The implementation of the Multiply class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "multiply.h"
#include "messages.h"
#include "logger.h"
#include "index.h"

/** The constructor.
 *
 * @param ANode The node of matrix A.
 * @param BNode The node of matrix B.
 * @param CNode The node of matrix C.
 */
MultiplyElement::MultiplyElement (int blocksize,
    CProxyElement_Node ANode,
    CProxyElement_Node BNode,
    CProxyElement_Node CNode)
{
  DEBUG("initializing multiply element\n");
  this->blocksize = blocksize;
  this->ANode = ANode;
  this->BNode = BNode;
  this->CNode = CNode;
}

/** The migration method.
 *
 * @param msg The migration message.
 */
MultiplyElement::MultiplyElement (CkMigrateMessage *msg)
{
  ABORT("migrating\n");
}

/** Multiply nodes.
 */
void MultiplyElement::multiply (CkCallback &done)
{
  DEBUG("(%d,%d,%d) multiply\n", thisIndex.x, thisIndex.y, thisIndex.z);

  NodeBlockMsg *ABlock = ANode.getBlock();
  NodeBlockMsg *BBlock = BNode.getBlock();
  NodeBlockMsg *CBlock = CNode.getBlock();

  for(int i = 0; i < blocksize; i++) {
    for(int j = 0; j < blocksize; j++) {
      for(int k = 0; k < blocksize; k++)
      {
        CBlock->block[BLOCK_INDEX(i, j, 0, 0, blocksize)] +=
          ABlock->block[BLOCK_INDEX(i, k, 0, 0, blocksize)]
          *BBlock->block[BLOCK_INDEX(k, j, 0, 0, blocksize)];
      }
    }
  }
  contribute(done);
}

/** The constructor.
 */
Multiply::Multiply ()
{
}

/** Multiply two Matrix objects.
 *
 * A full convolution (as a space filling curve) is constructed.
 *
 * @param A Matrix A.
 * @param B Matrix B.
 * @param C Matrix C.
 * @param cb The callback.
 */
void Multiply::multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C,
    CkCallback &cb)
{
  DEBUG("initializing multiply\n");

  MatrixInfoMsg *AInfo = A.info();
  MatrixInfoMsg *BInfo = B.info();
  MatrixInfoMsg *CInfo = C.info();

  if(AInfo->N != BInfo->N || AInfo->N != CInfo->N)
  {
    ABORT("matrix dimension mismatch\n");
  }

  if(AInfo->blocksize != BInfo->blocksize || AInfo->blocksize != CInfo->blocksize)
  {
    ABORT("blocksize mismatch\n");
  }

  CProxyElement_Node ANode = AInfo->tierNode(0, 0);

  convolution = CProxy_MultiplyElement::ckNew();
  for(int i = 0; i < (1 << CInfo->depth); i++) {
    for(int j = 0; j < (1 << CInfo->depth); j++) {
      for(int k = 0; k < (1 << CInfo->depth); k++)
      {
        DEBUG("inserting convolution at C(%d,%d) <- A(%d,%d) * B(%d,%d)\n",
            i, j, i, k, k, j);
        convolution(i, j, k).insert(CInfo->blocksize,
            AInfo->tierNode(i, k),
            BInfo->tierNode(k, j),
            CInfo->tierNode(i, j));
      }
    }
  }
  convolution.doneInserting();
  DEBUG("done initializing convolution\n");

  DEBUG("multiplying...\n");
  done = cb;
  CkCallback elementsDone(CkReductionTarget(Multiply, multiplyDone), thisProxy);
  convolution.multiply(elementsDone);
  DEBUG("done here\n");
}

/** Reduction target.
 */
void Multiply::multiplyDone ()
{
  DEBUG("done, sending back\n");
  done.send();
}

#include "multiply.def.h"
