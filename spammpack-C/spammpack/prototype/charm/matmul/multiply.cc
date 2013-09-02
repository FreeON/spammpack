/** @file
 *
 * The implementation of the Multiply class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"
#include "multiply.h"
#include "multiplyelement.h"
#include "messages.h"
#include "logger.h"
#include "utilities.h"
#include "index.h"

#include <sstream>

/** The constructor.
 *
 * @param A Matrix A.
 * @param B Matrix B.
 * @param C Matrix C.
 * @param blocksize The SpAMM blocksize.
 * @param depth The depth of the matrix trees.
 * @param ANodes The Node objects of A.
 * @param BNodes The Node objects of B.
 * @param CNodes The Node objects of C.
 */
Multiply::Multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C,
    int blocksize, int depth, CProxy_Node ANodes, CProxy_Node BNodes,
    CProxy_Node CNodes)
{
  this->A = A;
  this->B = B;
  this->C = C;

  this->depth = depth;

  convolution = new CProxy_MultiplyElement[depth+1];
  for(int tier = depth; tier >= 0; tier--)
  {
    int NTier = 1 << tier;

    unsigned long bytes = NTier*NTier*NTier*(sizeof(MultiplyElement)+blocksize*blocksize*sizeof(double));
    INFO("created %dx%dx%d convolution, %d MultiplyElements "
        "using %d bytes (%s)\n",
        NTier, NTier, NTier, NTier*NTier*NTier,
        bytes, humanReadableSize(bytes).c_str());

    convolution[tier] = CProxy_MultiplyElement::ckNew(blocksize, tier, depth,
        ANodes, BNodes, CNodes, NTier, NTier, NTier);

    if(tier < depth)
    {
      convolution[tier].setNextConvolution(convolution[tier+1]);
    }

    if(tier == depth)
    {
      PEMap = new int[NTier*NTier*NTier];
    }
  }
}

/** The destructor.
 */
Multiply::~Multiply (void)
{
  delete[] PEMap;
}

/** Multiply two Matrix objects.
 *
 * A full convolution (as a space filling curve) is constructed.
 *
 * @param tolerance The SpAMM tolerance.
 * @param cb The callback.
 */
void Multiply::multiply (double tolerance, CkCallback &cb)
{
  INFO("tolerance = %e\n", tolerance);

  for(int tier = 0; tier < depth; tier++)
  {
    /* Prune convolution to achieve reduced complexity in symbolic part of the
     * multiply. */
    INFO("pruning tier %d\n", tier+1);
    MatrixNodeMsg *ANodes = A.getNodes(tier+1);
    MatrixNodeMsg *BNodes = B.getNodes(tier+1);
    int NTier = 1 << (tier+1);
    convolution[tier].pruneProduct(tolerance, ANodes->nodes, BNodes->nodes,
        CkCallbackResumeThread());
    delete ANodes;
    delete BNodes;

#ifdef PRUNE_CONVOLUTION
    /* Mark the next tier as complete so that the load balancer can work. */
    INFO("calling doneInserting() on tier %d\n", tier+1);
    convolution[tier+1].doneInserting();

    /* Wait until things have settled down. */
    CkWaitQD();
#endif
  }

  /* Multiply. */
  convolution[depth].multiply(tolerance, CkCallbackResumeThread());

  INFO("storeBack\n");
  convolution[depth].storeBack(CkCallbackResumeThread());

  /* Update norms. */
  INFO("update norms\n");
  C.setNorm(CkCallbackResumeThread());

  INFO("done\n");
  cb.send();
}

/** Print the PEs all @link Node Nodes @endlink are on.
 *
 * @param cb The callback to signal once all @link Node Nodes @endlink have
 * printed.
 */
void Multiply::updatePEMap (CkCallback &cb)
{
  this->cb = cb;
  CkCallback done(CkReductionTarget(Multiply, donePEMap), thisProxy);
  convolution[depth].PEMap(done);
}

void Multiply::donePEMap (CkReductionMsg *msg)
{
  int NTier = 1 << depth;
  CkReduction::setElement *current = (CkReduction::setElement*) msg->getData();
  while(current != NULL)
  {
    int *intPtr = (int*) &current->data;
    DEBUG("data = { %d, %d, %d }\n", intPtr[0], intPtr[1], intPtr[2], intPtr[3]);
    PEMap[BLOCK_INDEX_3(intPtr[0], intPtr[1], intPtr[2], NTier)] = intPtr[3];
    current = current->next();
  }

#ifdef ATOMIC_PE_MAP
  std::ostringstream o;
  o << "PEMap for convolution:" << std::endl;
#else
  INFO("PEMap for convolution:\n");
#endif
  for(int i = 0; i < NTier; i++) {
    for(int j = 0; j < NTier; j++) {
      for(int k = 0; k < NTier; k++)
      {
#ifdef ATOMIC_PE_MAP
        o << "PEMap(" << i << "," << j << "," << k << ") = "
          << PEMap[BLOCK_INDEX_3(i, j, k, NTier)] << std::endl;;
#else
        CkPrintf("PEMap(%d,%d,%d) = %d\n", i, j, k, PEMap[BLOCK_INDEX_3(i, j, k, NTier)]);
#endif
      }
    }
  }
#ifdef ATOMIC_PE_MAP
  o << "end of PEMap for convolution" << std::endl;
  INFO(o.str().c_str());
#else
  CkPrintf("end of PEMap for convolution\n");
#endif

  cb.send();
}

#include "multiply.def.h"
