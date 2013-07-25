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
#include <string.h>

/** The constructor.
 *
 * @param ANode The node of matrix A.
 * @param BNode The node of matrix B.
 * @param CNode The node of matrix C.
 */
MultiplyElement::MultiplyElement (int blocksize, CProxy_Node A,
    CProxy_Node B, CProxy_Node C)
{
  DEBUG("initializing multiply element\n");
  this->blocksize = blocksize;
  this->ANode = A(thisIndex.x, thisIndex.z);
  this->BNode = B(thisIndex.z, thisIndex.y);
  this->CNode = C(thisIndex.x, thisIndex.y);
  CResult = NULL;
}

/** The destructor.
 */
MultiplyElement::~MultiplyElement ()
{
  delete[] CResult;
}

/** The migration method.
 *
 * @param msg The migration message.
 */
MultiplyElement::MultiplyElement (CkMigrateMessage *msg)
{
  ABORT("migrating\n");
}

/** The PUP method.
 */

void MultiplyElement::pup (PUP::er &p)
{
  p|index;
  p|blocksize;
  p|ANode;
  p|BNode;
  p|CNode;
  int CResultNull = (CResult == NULL);
  if(CResultNull)
  {
    CResult = NULL;
  }
  else
  {
    if(p.isUnpacking())
    {
      CResult = new double[blocksize*blocksize];
    }
    p|*CResult;
  }
}

/** Multiply nodes.
 */
void MultiplyElement::multiply (CkCallback &done)
{
  DEBUG("(%d,%d,%d) multiply\n", thisIndex.x, thisIndex.y, thisIndex.z);

  NodeBlockMsg *ABlock = ANode.getBlock();
  NodeBlockMsg *BBlock = BNode.getBlock();

  CResult = new double[blocksize*blocksize];
  memset(CResult, 0, sizeof(double)*blocksize*blocksize);

  for(int i = 0; i < blocksize; i++) {
    for(int j = 0; j < blocksize; j++) {
      for(int k = 0; k < blocksize; k++)
      {
        CResult[BLOCK_INDEX(i, j, 0, 0, blocksize)] +=
          ABlock->block[BLOCK_INDEX(i, k, 0, 0, blocksize)]
          *BBlock->block[BLOCK_INDEX(k, j, 0, 0, blocksize)];
      }
    }
  }
  contribute(done);
}

/** Push the C Nodes back into the C Matrix.
 */
void MultiplyElement::storeBack (CkCallback &done)
{
  CNode.add(blocksize, CResult);
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

  convolution = CProxy_MultiplyElement::ckNew();
  for(int i = 0; i < (1 << CInfo->depth); i++) {
    for(int j = 0; j < (1 << CInfo->depth); j++)
    {
      CProxyElement_Node CNode = CInfo->tierNode(i, j);
      for(int k = 0; k < (1 << CInfo->depth); k++)
      {
        CProxyElement_Node ANode = AInfo->tierNode(i, k);
        CProxyElement_Node BNode = BInfo->tierNode(k, j);

        DEBUG("inserting convolution at C(%d,%d) <- A(%d,%d) * B(%d,%d)\n",
            i, j, i, k, k, j);

        convolution(i, j, k).insert(CInfo->blocksize, AInfo->tierNode,
            BInfo->tierNode, CInfo->tierNode);

        DEBUG("inserted element\n");
      }
    }
  }
  convolution.doneInserting();
  DEBUG("done initializing convolution\n");

  DEBUG("multiplying...\n");
  convolution.multiply(CkCallbackResumeThread());
  convolution.storeBack(CkCallbackResumeThread());
  cb.send();
}

#include "multiply.def.h"
