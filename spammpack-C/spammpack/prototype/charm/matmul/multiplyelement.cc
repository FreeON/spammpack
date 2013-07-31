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
MultiplyElement::MultiplyElement (int blocksize, int tier, int depth,
    CProxy_Node A, CProxy_Node B, CProxy_Node C)
{
  DEBUG("tier %d ME(%d,%d,%d) constructor\n", tier, thisIndex.x, thisIndex.y,
      thisIndex.z);

  this->blocksize = blocksize;
  this->tier = tier;
  this->depth = depth;
  this->A = A;
  this->B = B;
  this->C = C;

  if(tier < depth)
  {
    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++) {
        for(int k = 0; k < 2; k++)
        {
          convolutionExists[i][j][k] = true;
        }
      }
    }
  }

  CResult = NULL;
  numberCalls = 0;
}

/** The migration constructor.
 *
 * @param msg The migration message.
 */
MultiplyElement::MultiplyElement (CkMigrateMessage *msg)
{
  DEBUG("tier %d ME(%d,%d,%d) migration constructor\n", tier, thisIndex.x,
      thisIndex.y, thisIndex.z);
}

/** The destructor.
 */
MultiplyElement::~MultiplyElement ()
{
  DEBUG("tier %d ME(%d,%d,%d) destructor\n", tier, thisIndex.x, thisIndex.y,
      thisIndex.z);
  delete[] CResult;

  if(tier < depth)
  {
    /* Recursively prune convolution elements. */
    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++) {
        for(int k = 0; k < 2; k++)
        {
          if(convolutionExists[i][j][k])
          {
            nextConvolution((thisIndex.x << 1)+i, (thisIndex.y << 1)+j,
                (thisIndex.z << 1)+k).ckDestroy();
          }
        }
      }
    }
  }
}

/** The the next convolution array.
 *
 * @param nextConvolution The next convolution array.
 */
void MultiplyElement::setNextTier (CProxy_MultiplyElement nextConvolution,
    CProxy_Node nextA, CProxy_Node nextB, CkCallback &cb)
{
  this->nextConvolution = nextConvolution;
  this->nextA = nextA;
  this->nextB = nextB;
  contribute(cb);
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
  p|depth;
  p|tier;
  p|A;
  p|B;
  p|C;

  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++) {
      for(int k = 0; k < 2; k++)
      {
        p|convolutionExists[i][j][k];
      }
    }
  }

  if(tier < depth)
  {
    p|nextConvolution;
    p|nextA;
    p|nextB;
  }

  int numberElements = (CResult == NULL ? 0 : blocksize*blocksize);
  p|numberElements;

  if(p.isUnpacking())
  {
    DEBUG("tier %d ME(%d,%d,%d) pup: unpacking %d elements\n",
        tier, thisIndex.x, thisIndex.y, thisIndex.z, numberElements);
  }
  else
  {
    if(p.isSizing())
    {
      DEBUG("tier %d ME(%d,%d,%d) pup: sizing %d elements\n",
          tier, thisIndex.x, thisIndex.y, thisIndex.z, numberElements);
    }
    else
    {
      DEBUG("tier %d ME(%d,%d,%d) pup: packing %d elements\n",
          tier, thisIndex.x, thisIndex.y, thisIndex.z, numberElements);
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
void MultiplyElement::multiply (double tolerance, CkCallback &cb)
{
  DEBUG("tier %d ME(%d,%d,%d) multiply\n", tier, thisIndex.x, thisIndex.y,
      thisIndex.z);

  if(numberCalls > 0)
  {
    ABORT("tier %d ME(%d,%d,%d) this MultiplyElement has been called before\n",
        tier, thisIndex.x, thisIndex.y, thisIndex.z);
  }
  numberCalls++;

  if(tier == depth)
  {
    NodeInfoMsg *AInfo = A(thisIndex.x, thisIndex.z).info();
    NodeInfoMsg *BInfo = B(thisIndex.z, thisIndex.y).info();

    if(AInfo->norm*BInfo->norm > tolerance)
    {
      if(CResult != NULL)
      {
        ABORT("tier %d ME(%d,%d,%d) CResult is not NULL\n", tier, thisIndex.x,
            thisIndex.y, thisIndex.z);
      }

      DEBUG("tier %d ME(%d,%d,%d) multiplying blocks\n", tier, thisIndex.x,
          thisIndex.y, thisIndex.z);

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
    }

    else
    {
      ABORT("tier %d ME(%d,%d,%d) skipping block produce\n", tier,
          thisIndex.x, thisIndex.y, thisIndex.z);
    }

    delete AInfo;
    delete BInfo;
  }

  else
  {
    /* Get information on A and B matrices. */
    NodeInfoMsg *AInfo[2][2];
    NodeInfoMsg *BInfo[2][2];

    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++)
      {
        AInfo[i][j] = A(thisIndex.x+i, thisIndex.z+j).info();
        BInfo[i][j] = B(thisIndex.z+i, thisIndex.y+j).info();
      }
    }

    /* Check what products are necessary one tier down. */
    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++) {
        for(int k = 0; k < 2; k++)
        {
          if(AInfo[i][k]->norm*BInfo[k][j]->norm > tolerance)
          {
            if(!convolutionExists[i][j][k])
            {
              DEBUG("tier %d ME(%d,%d,%d) adding product C(%d,%d) <- A(%d,%d)*B(%d,%d)\n",
                  tier, thisIndex.x, thisIndex.y, thisIndex.z, i, j, i, k, k, j);
              nextConvolution((thisIndex.x << 1)+i, (thisIndex.y << 1)+j,
                  (thisIndex.z << 1)+k).insert(blocksize, tier+1, depth, A, B, C);
              convolutionExists[i][j][k] = true;
            }

            else
            {
              DEBUG("tier %d ME(%d,%d,%d) keeping product C(%d,%d) <- A(%d,%d)*B(%d,%d)\n",
                  tier, thisIndex.x, thisIndex.y, thisIndex.z, i, j, i, k, k, j);
            }
          }

          else
          {
            if(convolutionExists[i][j][k])
            {
              DEBUG("tier %d ME(%d,%d,%d) dropping product C(%d,%d) <- A(%d,%d)*B(%d,%d)\n",
                  tier, thisIndex.x, thisIndex.y, thisIndex.z, i, j, i, k, k, j);
              nextConvolution((thisIndex.x << 1)+i, (thisIndex.y << 1)+j,
                  (thisIndex.z << 1)+k).ckDestroy();
              convolutionExists[i][j][k] = false;
            }

            else
            {
              DEBUG("tier %d ME(%d,%d,%d) product already does not exist C(%d,%d) <- A(%d,%d)*B(%d,%d)\n",
                  tier, thisIndex.x, thisIndex.y, thisIndex.z, i, j, i, k, k, j);
            }
          }
        }
      }
    }

    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++)
      {
        delete AInfo[i][j];
        delete BInfo[i][j];
      }
    }
  }

  DEBUG("tier %d ME(%d,%d,%d) contribute\n", tier, thisIndex.x, thisIndex.y,
      thisIndex.z);
  contribute(cb);
}

/** Push the C submatrices back into the C Matrix.
 *
 * @param cb The callback.
 */
void MultiplyElement::storeBack (CkCallback &cb)
{
#ifdef DEBUG_OUTPUT
  DEBUG("ME(%d,%d,%d) storing back\n", thisIndex.x, thisIndex.y, thisIndex.z);
  printDense(blocksize, CResult);
#endif
  C(thisIndex.x, thisIndex.y).add(blocksize, CResult);
  contribute(cb);
}

#include "multiplyelement.def.h"
