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
  this->A = A;
  this->B = B;
  this->C = C;
  CResult = NULL;
  numberCalls = 0;
}

/** The migration constructor.
 *
 * @param msg The migration message.
 */
MultiplyElement::MultiplyElement (CkMigrateMessage *msg)
{
  INFO("ME(%d,%d,%d) migration constructor\n", thisIndex.x, thisIndex.y, thisIndex.z);
}

/** The destructor.
 */
MultiplyElement::~MultiplyElement ()
{
  INFO("ME(%d,%d,%d) destructor\n", thisIndex.x, thisIndex.y, thisIndex.z);
  delete[] CResult;
  INFO("ME(%d,%d,%d) destructor done\n", thisIndex.x, thisIndex.y, thisIndex.z);
}

/** The PUP method.
 */

void MultiplyElement::pup (PUP::er &p)
{
  p|index;
  p|blocksize;
  p|A;
  p|B;
  p|C;

  int numberElements = (CResult == NULL ? 0 : blocksize*blocksize);
  p|numberElements;

  if(p.isUnpacking())
  {
    INFO("pup: ME(%d,%d,%d) unpacking %d elements\n",
        thisIndex.x, thisIndex.y, thisIndex.z, numberElements);
  }
  else
  {
    if(p.isSizing())
    {
      INFO("pup: ME(%d,%d,%d) sizing %d elements\n",
          thisIndex.x, thisIndex.y, thisIndex.z, numberElements);
    }
    else
    {
      INFO("pup: ME(%d,%d,%d) packing %d elements\n",
          thisIndex.x, thisIndex.y, thisIndex.z, numberElements);
    }
  }

  if(numberElements > 0)
  {
    if(p.isUnpacking())
    {
      CResult = new double[numberElements];
    }
    PUParray(p, CResult, numberElements);
  }
  else
  {
    if(p.isUnpacking()) { CResult = NULL; }
  }

  p|numberCalls;
}

/** Multiply nodes.
 */
void MultiplyElement::multiply ()
{
  INFO("ME(%d,%d,%d) multiply\n", thisIndex.x, thisIndex.y, thisIndex.z);

  if(numberCalls > 0)
  {
    ABORT("ME(%d,%d,%d) this MultiplyElement has been called before\n",
        thisIndex.x, thisIndex.y, thisIndex.z);
  }
  numberCalls++;

  if(CResult != NULL)
  {
    ABORT("ME(%d,%d,%d) CResult is not NULL\n", thisIndex.x, thisIndex.y, thisIndex.z);
  }

  CResult = new double[blocksize*blocksize];
  memset(CResult, 0, sizeof(double)*blocksize*blocksize);

  NodeBlockMsg *ABlock = A(thisIndex.x, thisIndex.z).getBlock();
  NodeBlockMsg *BBlock = B(thisIndex.z, thisIndex.y).getBlock();

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
  contribute();
  migrateMe(0);
}

/** Push the C Nodes back into the C Matrix.
 */
void MultiplyElement::storeBack ()
{
  INFO("ME(%d,%d,%d) storing back\n", thisIndex.x, thisIndex.y, thisIndex.z);
  C(thisIndex.x, thisIndex.y).add(blocksize, CResult);
  contribute();
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
  this->cb = cb;
  convolution.ckSetReductionClient(new CkCallback(CkReductionTarget(Multiply, multiplyDone), thisProxy));
  convolution.multiply();
}

/** Reduction target for MultiplyElement::multiply().
 */
void Multiply::multiplyDone ()
{
  INFO("multiplyDone\n");
  convolution.ckSetReductionClient(new CkCallback(CkReductionTarget(Multiply, storeBackDone), thisProxy));
  convolution.storeBack();
}

void Multiply::storeBackDone ()
{
  INFO("storeBackDone\n");
  cb.send();
}

#include "multiply.def.h"
