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
 *
 * A full convolution (as a space filling curve) is constructed.
 *
 * @param A Matrix A.
 * @param B Matrix B.
 * @param C Matrix C.
 */
Multiply::Multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C)
{
  DEBUG("initializing multiply\n");

  MatrixInfoMsg *AInfo = A.info();
  DEBUG("here\n");
  MatrixInfoMsg *BInfo = B.info();
  DEBUG("here\n");
  MatrixInfoMsg *CInfo = C.info();
  DEBUG("here\n");

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
    for(int j = 0; j < (1 << CInfo->depth); j++) {
      for(int k = 0; k < (1 << CInfo->depth); k++)
      {
        DEBUG("inserting convolution at (%d,%d) <- (%d,%d) * (%d,%d)\n",
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
}

/** Multiply two Matrix objects.
 *
 * @param cb The callback.
 */
void Multiply::multiply (CkCallback &cb)
{
  DEBUG("multiplying...\n");
  done = cb;
  CkCallback elementsDone(CkReductionTarget(Multiply, multiplyDone), thisProxy);
  convolution.multiply(elementsDone);
}

/** Reduction target.
 */
void Multiply::multiplyDone ()
{
  done.send();
}

//void Node::multiply (int index, CProxy_Node A, CProxy_Node B, CkCallback &cb)
//{
//  DEBUG("(%d) starting\n", index);
//
//  callbackSet[index] = true;
//  this->cb[index] = cb;
//
//  if(tier == depth)
//  {
//    DEBUG("(%d) reached depth\n", index);
//    NodeBlockMsg *ABlock = A.getBlock();
//    NodeBlockMsg *BBlock = B.getBlock();
//
//    for(int i = 0; i < blocksize; i++) {
//      for(int j = 0; j < blocksize; j++) {
//        for(int k = 0; k < blocksize; k++)
//        {
//          block[BLOCK_INDEX(i, j, 0, 0, blocksize)] +=
//            ABlock->block[BLOCK_INDEX(i, k, 0, 0, blocksize)]
//            *BBlock->block[BLOCK_INDEX(k, j, 0, 0, blocksize)];
//        }
//      }
//    }
//
//    cb.send(new IntMsg(index));
//  }
//
//  else
//  {
//    NodeInfoMsg *AInfo = A.info();
//    NodeInfoMsg *BInfo = B.info();
//
//    for(int i = 0; i < 2; i++) {
//      for(int j = 0; j < 2; j++)
//      {
//        if(childNull[CHILD_INDEX(i, j)])
//        {
//          ABORT("create new child\n");
//        }
//
//        for(int k = 0; k < 2; k++)
//        {
//          int childIndex = (index << 3) | (i << 2) | (j << 1) | k;
//          childWorking[index].push_back(childIndex);
//          CkCallback thisCB(CkIndex_Node::multiplyDone(NULL), thisProxy);
//          DEBUG("(%d) descending C(%d,%d) <- A(%d,%d)*B(%d,%d)\n", index,
//              i, j, i, k, k, j);
//          child[CHILD_INDEX(i, j)].multiply(childIndex,
//              AInfo->child[CHILD_INDEX(i, k)],
//              BInfo->child[CHILD_INDEX(k, j)],
//              thisCB);
//        }
//      }
//    }
//  }
//}
//
//void Node::multiplyDone (IntMsg *index)
//{
//  int thisIndex = index->i >> 3;
//  DEBUG("(%d) received index %d\n", thisIndex, index->i);
//  for(std::list<int>::iterator i = childWorking[thisIndex].begin();
//      i != childWorking[thisIndex].end();
//      i++)
//  {
//    if(*i == index->i)
//    {
//      childWorking[thisIndex].erase(i);
//      break;
//    }
//  }
//
//  if(childWorking[thisIndex].size() > 0)
//  {
//    return;
//  }
//
//  DEBUG("(%d) sending to callback\n", thisIndex);
//
//  childWorking.erase(thisIndex);
//  callbackSet[index->i] = false;
//  CkCallback tempCb = cb[thisIndex];
//  cb.erase(thisIndex);
//  tempCb.send(new IntMsg(thisIndex));
//}

#include "multiply.def.h"
