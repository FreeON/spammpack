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
          nextConvolutionExists[i][j][k] = true;
        }
      }
    }
  }

  CResult = NULL;
  numberCalls = 0;
  wasMigrated = false;
}

/** The migration constructor.
 *
 * @param msg The migration message.
 */
MultiplyElement::MultiplyElement (CkMigrateMessage *msg)
{
  INFO("ME(%d,%d,%d) migration constructor\n", thisIndex.x, thisIndex.y,
      thisIndex.z);

  /* Reset the migration flag on the new PE. */
  wasMigrated = false;
}

/** The destructor.
 */
MultiplyElement::~MultiplyElement ()
{
  DEBUG("tier %d ME(%d,%d,%d) destructor\n", tier, thisIndex.x, thisIndex.y,
      thisIndex.z);

  delete[] CResult;

  if(!wasMigrated && (tier < depth))
  {
    DEBUG("tier %d ME(%d,%d,%d) destructor, pruning lower tier elements\n",
        tier, thisIndex.x, thisIndex.y, thisIndex.z);

    /* Recursively prune convolution elements. */
    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++) {
        for(int k = 0; k < 2; k++)
        {
          if(nextConvolutionExists[i][j][k])
          {
            DEBUG("tier %d ME(%d,%d,%d) destructor, pruning tier %d ME(%d,%d,%d)\n",
                tier, thisIndex.x, thisIndex.y, thisIndex.z,
                tier+1,
                (thisIndex.x << 1)+i,
                (thisIndex.y << 1)+j,
                (thisIndex.z << 1)+k);

            nextConvolution((thisIndex.x << 1)+i, (thisIndex.y << 1)+j,
                (thisIndex.z << 1)+k).ckDestroy();

            DEBUG("tier %d ME(%d,%d,%d) destructor, pruning tier %d ME(%d,%d,%d), here 2\n",
                tier, thisIndex.x, thisIndex.y, thisIndex.z,
                tier+1,
                (thisIndex.x << 1)+i,
                (thisIndex.y << 1)+j,
                (thisIndex.z << 1)+k);

            nextConvolutionExists[i][j][k] = false;

            DEBUG("tier %d ME(%d,%d,%d) destructor, pruning tier %d ME(%d,%d,%d), here 3\n",
                tier, thisIndex.x, thisIndex.y, thisIndex.z,
                tier+1,
                (thisIndex.x << 1)+i,
                (thisIndex.y << 1)+j,
                (thisIndex.z << 1)+k);
          }
        }
      }
    }
  }

  DEBUG("tier %d ME(%d,%d,%d) destructor done\n", tier, thisIndex.x,
      thisIndex.y, thisIndex.z);
}

/** The the next convolution array.
 *
 * @param nextConvolution The next convolution array.
 */
void MultiplyElement::setNextTier (CProxy_MultiplyElement nextConvolution,
    CProxy_Node nextA, CProxy_Node nextB, CkCallback &cb)
{
  DEBUG("tier %d ME(%d,%d,%d) setting nextTier\n", tier, thisIndex.x, thisIndex.y, thisIndex.z);
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
        p|nextConvolutionExists[i][j][k];
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
    INFO("tier %d ME(%d,%d,%d) pup: unpacking %d elements\n",
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
      INFO("tier %d ME(%d,%d,%d) pup: packing %d elements\n",
          tier, thisIndex.x, thisIndex.y, thisIndex.z, numberElements);

      /* Set the wasMigrated flag to indicate that this instance is going to
       * get destroyed because of a migration, and not because of pruning. */
      wasMigrated = true;

      print("packing");
    }
  }

  if(numberElements > 0)
  {
    if(p.isUnpacking())
    {
      CResult = new double[numberElements];
    }
    PUParray(p, CResult, numberElements);

    print("unpacking");
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

      INFO("tier %d ME(%d,%d,%d) sleeping\n", tier, thisIndex.x,
          thisIndex.y, thisIndex.z);
      sleep(30);
    }

    else
    {
      ABORT("tier %d ME(%d,%d,%d) skipping block product\n", tier,
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
        AInfo[i][j] = A((thisIndex.x << 1)+i, (thisIndex.z << 1)+j).info();
        BInfo[i][j] = B((thisIndex.z << 1)+i, (thisIndex.y << 1)+j).info();
      }
    }

    /* Check what products are necessary one tier down. */
    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++) {
        for(int k = 0; k < 2; k++)
        {
          if(AInfo[i][k]->norm*BInfo[k][j]->norm > tolerance)
          {
            if(!nextConvolutionExists[i][j][k])
            {
              DEBUG("tier %d ME(%d,%d,%d) adding product ME(%d,%d,%d): C(%d,%d) <- A(%d,%d)*B(%d,%d)\n",
                  tier, thisIndex.x, thisIndex.y, thisIndex.z,
                  (thisIndex.x << 1)+i, (thisIndex.y << 1)+j, (thisIndex.z << 1)+k,
                  (thisIndex.x << 1)+i, (thisIndex.y << 1)+j,
                  (thisIndex.x << 1)+i, (thisIndex.z << 1)+k,
                  (thisIndex.z << 1)+k, (thisIndex.y << 1)+j);
              nextConvolution((thisIndex.x << 1)+i, (thisIndex.y << 1)+j,
                  (thisIndex.z << 1)+k).insert(blocksize, tier+1, depth, A, B, C);
              nextConvolutionExists[i][j][k] = true;
            }

            else
            {
              DEBUG("tier %d ME(%d,%d,%d) keeping product ME(%d,%d,%d): C(%d,%d) <- A(%d,%d)*B(%d,%d)\n",
                  tier, thisIndex.x, thisIndex.y, thisIndex.z,
                  (thisIndex.x << 1)+i, (thisIndex.y << 1)+j, (thisIndex.z << 1)+k,
                  (thisIndex.x << 1)+i, (thisIndex.y << 1)+j,
                  (thisIndex.x << 1)+i, (thisIndex.z << 1)+k,
                  (thisIndex.z << 1)+k, (thisIndex.y << 1)+j);
            }
          }

          else
          {
            if(nextConvolutionExists[i][j][k])
            {
              DEBUG("tier %d ME(%d,%d,%d) dropping product ME(%d,%d,%d): C(%d,%d) <- A(%d,%d)*B(%d,%d)\n",
                  tier, thisIndex.x, thisIndex.y, thisIndex.z,
                  (thisIndex.x << 1)+i, (thisIndex.y << 1)+j, (thisIndex.z << 1)+k,
                  (thisIndex.x << 1)+i, (thisIndex.y << 1)+j,
                  (thisIndex.x << 1)+i, (thisIndex.z << 1)+k,
                  (thisIndex.z << 1)+k, (thisIndex.y << 1)+j);

              nextConvolutionExists[i][j][k] = false;
              nextConvolution((thisIndex.x << 1)+i, (thisIndex.y << 1)+j,
                  (thisIndex.z << 1)+k).ckDestroy();
            }

            else
            {
              DEBUG("tier %d ME(%d,%d,%d) product already does not exist ME(%d,%d,%d): C(%d,%d) <- A(%d,%d)*B(%d,%d)\n",
                  tier, thisIndex.x, thisIndex.y, thisIndex.z,
                  (thisIndex.x << 1)+i, (thisIndex.y << 1)+j, (thisIndex.z << 1)+k,
                  (thisIndex.x << 1)+i, (thisIndex.y << 1)+j,
                  (thisIndex.x << 1)+i, (thisIndex.z << 1)+k,
                  (thisIndex.z << 1)+k, (thisIndex.y << 1)+j);
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

#ifdef FORCED_MIGRATION
  if(CkNumPes() > 1)
  {
    INFO("tier %d ME(%d,%d,%d) requesting migration to PE %d\n",
        tier, thisIndex.x, thisIndex.y, thisIndex.z, (CkMyPe()+1)%CkNumPes());
    migrateMe((CkMyPe()+1)%CkNumPes());
  }
#endif
}

/** Push the C submatrices back into the C Matrix.
 *
 * @param cb The callback.
 */
void MultiplyElement::storeBack (CkCallback &cb)
{
#ifdef DEBUG_OUTPUT
  DEBUG("tier %d ME(%d,%d,%d) storing back\n", tier, thisIndex.x, thisIndex.y, thisIndex.z);
  printDense(blocksize, CResult);
#endif
  C(thisIndex.x, thisIndex.y).add(blocksize, CResult);
  contribute(cb);
}

/** Print a MultiplyElement.
 */
void MultiplyElement::print (std::string tag)
{
  INFO("tier %d ME(%d,%d,%d) [%s] CResult = %p\n", tier, thisIndex.x, thisIndex.y,
      thisIndex.z, tag.c_str(), CResult);
}

#include "multiplyelement.def.h"
