/** @file
 *
 * The implementation of the MultiplyElement class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "multiplyelement.h"
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
  DEBUG("ME(%d,%d,%d) constructor\n", thisIndex.x, thisIndex.y, thisIndex.z);
  this->blocksize = blocksize;
  this->A = A;
  this->B = B;
  this->C = C;
  CResult = NULL;
  numberCalls = 0;
  DEBUG("ME(%d,%d,%d) constructor done\n", thisIndex.x, thisIndex.y, thisIndex.z);
}

/** The migration constructor.
 *
 * @param msg The migration message.
 */
MultiplyElement::MultiplyElement (CkMigrateMessage *msg)
{
  DEBUG("ME(%d,%d,%d) migration constructor\n", thisIndex.x, thisIndex.y, thisIndex.z);
}

/** The destructor.
 */
MultiplyElement::~MultiplyElement ()
{
  DEBUG("ME(%d,%d,%d) destructor\n", thisIndex.x, thisIndex.y, thisIndex.z);
  delete[] CResult;
  DEBUG("ME(%d,%d,%d) destructor done\n", thisIndex.x, thisIndex.y, thisIndex.z);
}

/** The PUP method.
 *
 * @param p The object.
 */
void MultiplyElement::pup (PUP::er &p)
{
  CBase_MultiplyElement::pup(p);
  p|index;
  p|blocksize;
  p|numberCalls;
  p|A;
  p|B;
  p|C;

  int numberElements = (CResult == NULL ? 0 : blocksize*blocksize);
  p|numberElements;

  if(p.isUnpacking())
  {
    DEBUG("ME(%d,%d,%d) pup: unpacking %d elements\n",
        thisIndex.x, thisIndex.y, thisIndex.z, numberElements);
  }
  else
  {
    if(p.isSizing())
    {
      DEBUG("ME(%d,%d,%d) pup: sizing %d elements\n",
          thisIndex.x, thisIndex.y, thisIndex.z, numberElements);
    }
    else
    {
      DEBUG("ME(%d,%d,%d) pup: packing %d elements\n",
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
}

/** Multiply nodes.
 *
 * @param cb The callback.
 */
void MultiplyElement::multiply (CkCallback &cb)
{
  DEBUG("ME(%d,%d,%d) multiply\n", thisIndex.x, thisIndex.y, thisIndex.z);

  if(numberCalls > 0)
  {
    ABORT("ME(%d,%d,%d) this MultiplyElement has been called before\n",
        thisIndex.x, thisIndex.y, thisIndex.z);
  }
  numberCalls++;

  DEBUG("ME(%d,%d,%d) here\n", thisIndex.x, thisIndex.y, thisIndex.z);

  if(CResult != NULL)
  {
    ABORT("ME(%d,%d,%d) CResult is not NULL\n", thisIndex.x, thisIndex.y, thisIndex.z);
  }

  CResult = new double[blocksize*blocksize];
  memset(CResult, 0, sizeof(double)*blocksize*blocksize);

  DEBUG("ME(%d,%d,%d) here\n", thisIndex.x, thisIndex.y, thisIndex.z);

  NodeBlockMsg *ABlock = A(thisIndex.x, thisIndex.z).getBlock();
  NodeBlockMsg *BBlock = B(thisIndex.z, thisIndex.y).getBlock();

  DEBUG("ME(%d,%d,%d) here\n", thisIndex.x, thisIndex.y, thisIndex.z);

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
  DEBUG("ME(%d,%d,%d) contribute\n", thisIndex.x, thisIndex.y, thisIndex.z);
  contribute(cb);
  DEBUG("ME(%d,%d,%d) migrate request\n", thisIndex.x, thisIndex.y, thisIndex.z);
  migrateMe(0);
}

/** Push the C submatrices back into the C Matrix.
 *
 * @param cb The callback.
 */
void MultiplyElement::storeBack (CkCallback &cb)
{
  DEBUG("ME(%d,%d,%d) storing back\n", thisIndex.x, thisIndex.y, thisIndex.z);
  C(thisIndex.x, thisIndex.y).add(blocksize, CResult);
  contribute(cb);
}

#include "multiplyelement.def.h"
