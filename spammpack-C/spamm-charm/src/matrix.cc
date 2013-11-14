/** @file
 *
 * The implementation of the Matrix class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "matrix.h"
#include "messages.h"
#include "logger.h"
#include "utilities.h"
#include "index.h"

#include <assert.h>

/** The constructor.
 *
 * @param initialPE The PE to place the Node chares.
 * @param alignPEs Align PEs in the diagonal matrix case.
 * @param N The matrix size.
 * @param blocksize The SpAMM blocksize.
 * @param nameLength The strlen of the name.
 * @param name The matrix name.
 */
Matrix::Matrix (int initialPE, bool alignPEs, int N, int blocksize,
    int nameLength, char *name)
{
  this->name = strdup(name);
  this->N = N;
  this->blocksize = blocksize;

  /* Calculate tree depth. */
  depth = -1;
  for(int i = N/blocksize; i > 0; i >>= 1)
  {
    depth++;
  }
  if(blocksize*(1 << depth) < N) depth++;
  NPadded = blocksize*(1 << depth);

  nodes = new CProxy_Node[depth+1];
  memset(nodes, 0, sizeof(CProxy_Node)*(depth+1));

  for(int tier = 0; tier <= depth; tier++)
  {
    int NTier = 1 << tier;

    unsigned long bytes = NTier*NTier*(sizeof(Node)
        +blocksize*blocksize*sizeof(double));
    INFO("name = %s, N = %d, blocksize = %d, tier = %d, depth = %d, "
        "NPadded = %d, NTier = %d, creating %d Nodes using %d bytes (%s)\n",
        this->name, N, blocksize, tier, depth, NPadded,
        NTier, NTier*NTier, bytes, humanReadableSize(bytes).c_str());

    nodes[tier] = CProxy_Node::ckNew();
    for(int i = 0; i < NTier; i++) {
      for(int j = 0; j < NTier; j++)
      {
        if(alignPEs && i == j)
        {
          initialPE = i%CkNumPes();
        }
        nodes[tier](i, j).insert(N, depth, blocksize, tier, initialPE);
      }
    }
    nodes[tier].doneInserting();

    if(tier == depth)
    {
      PEMap = new int[NTier*NTier];
      PEMap_norm = new double[NTier*NTier];
    }
  }
  DEBUG("done\n");
}

/** The destructor.
 */
Matrix::~Matrix (void)
{
  delete[] PEMap;
  delete[] PEMap_norm;
}

/** Initialize the matrix and its Nodes.
 *
 * @param cb The callback to notify once done.
 */
void Matrix::init (CkCallback &cb)
{
  INFO("initializing\n");
  for(int tier = 0; tier < depth+1; tier++)
  {
    nodes[tier].init(CkCallbackResumeThread());
  }
  cb.send();
}

/** Get some basic information on the matrix.
 *
 * @return The matrix information.
 */
MatrixInfoMsg * Matrix::info (void)
{
  return new MatrixInfoMsg (N, blocksize, depth, NPadded);
}

/** Convert a Matrix to a dense matrix.
 *
 * @return The dense matrix.
 */
DenseMatrixMsg * Matrix::toDense (void)
{
  DenseMatrixMsg *A = new (N*N) DenseMatrixMsg(N, N);

  int NTier = 1 << depth;

  for(int i = 0; i < NTier; i++) {
    for(int j = 0; j < NTier; j++)
    {
      DenseMatrixMsg *block = nodes[depth](i, j).toDense();

      for(int l = i*blocksize; l < (i+1)*blocksize && l < N; l++) {
        for(int m = j*blocksize; m < (j+1)*blocksize && m < N; m++)
        {
          A->A[BLOCK_INDEX(l, m, 0, 0, N)] = block->A[BLOCK_INDEX(l, m,
              i*blocksize, j*blocksize, blocksize)];
        }
      }

      delete block;
    }
  }

  A->M = N;
  A->N = N;

  return A;
}

/** Get the Node array on a particular tier.
 *
 * @return The Node array of the matrix.
 */
MatrixNodeMsg * Matrix::getNodes (void)
{
  DEBUG("getting nodes of matrix %s\n", name);
  MatrixNodeMsg *msg = new (depth+1) MatrixNodeMsg(depth+1);
  for(int tier = 0; tier < depth+1; tier++)
  {
    DEBUG("setting tier %d, isDelegated %d\n", tier, nodes[tier].ckIsDelegated());
    msg->nodes[tier] = nodes[tier];
  }
  return msg;
}

/** Print the PEs all @link Node Nodes @endlink are on.
 *
 * @param cb The callback to signal once all @link Node Nodes @endlink have
 * printed.
 */
void Matrix::updatePEMap (CkCallback &cb)
{
  this->cb = cb;
  CkCallback done(CkReductionTarget(Matrix, donePEMap), thisProxy);
  nodes[depth].PEMap(done);
}

/** The reduction target for Matrix::updatePEMap.
 *
 * @param msg The reduction message.
 */
void Matrix::donePEMap (CkReductionMsg *msg)
{
  int NTier = 1 << depth;

  CkReduction::setElement *current = (CkReduction::setElement*) msg->getData();
  while(current != NULL)
  {
    DEBUG("dataSize %d sizeof() %d\n", current->dataSize, sizeof(struct PEMap_Node_t));
    assert(current->dataSize == sizeof(struct PEMap_Node_t));
    struct PEMap_Node_t *result = (struct PEMap_Node_t*) &current->data;
    DEBUG("data = { %d, %d, %d }\n", result->index[0], result->index[1], result->PE);
    PEMap[BLOCK_INDEX(result->index[0], result->index[1], 0, 0, NTier)] = result->PE;
    PEMap_norm[BLOCK_INDEX(result->index[0], result->index[1], 0, 0, NTier)] = result->norm;
    current = current->next();
  }

  cb.send();
}

/** Set a matrix using a dense array.
 *
 * @param N The matrix size.
 * @param A The dense matrix.
 * @param cb The callback to signal once done.
 */
void Matrix::set (int N, double *A, CkCallback &cb)
{
  assert(this->N == N);

  DEBUG("setting %dx%d matrix %s\n", N, N, name);

  /* Set the A matrix. */
  double *block = new double[blocksize*blocksize];

  for(int i = 0; i < NPadded/blocksize; i++) {
    for(int j = 0; j < NPadded/blocksize; j++)
    {
      memset(block, 0, sizeof(double)*blocksize*blocksize);

      for(int l = i*blocksize; l < (i+1)*blocksize && l < N; l++) {
        for(int m = j*blocksize; m < (j+1)*blocksize && m < N; m++)
        {
          block[BLOCK_INDEX(l, m, i*blocksize, j*blocksize, blocksize)] =
            A[BLOCK_INDEX(l, m, 0, 0, N)];
        }
      }

#ifdef PRINT_MATRICES
      printDense(blocksize, block, "block(%i,%i)", i, j);
#endif

      DEBUG("calling set on Node(%d,%d)\n", i, j);
      nodes[depth](i, j).set(blocksize, block, CkCallbackResumeThread());
    }
  }

  delete[] block;

  /* Update norms. */
  thisProxy.setNorm(CkCallbackResumeThread());

  DEBUG("done setting matrix\n");
  cb.send();
}

/** Update the norms based on the norm information of the leaf @link Node
 * nodes @endlink.
 *
 * @param cb The callback to send back to.
 */
void Matrix::setNorm (CkCallback &cb)
{
  DEBUG("setting norm\n");

  /* Update the norms on the upper tiers. */
  for(int tier = depth-1; tier >= 0; tier--)
  {
    nodes[tier].setNorm(nodes[tier+1], CkCallbackResumeThread());
  }
  cb.send();
}

/** Return the PEMap. Call udpatePEMap() first.
 *
 * @return a PEMapMsg with the PEMap.
 */
PEMapMsg * Matrix::getPEMap (void)
{
  int NTier = 1 << depth;
  PEMapMsg *msg = new (NTier*NTier, NTier*NTier) PEMapMsg();
  memcpy(msg->PEMap, PEMap, NTier*NTier*sizeof(int));
  memcpy(msg->PEMap_norm, PEMap_norm, NTier*NTier*sizeof(double));
  return msg;
}

/** Update the trace on this matrix.
 *
 * @param cb The callback to send to when done.
 */
void Matrix::updateTrace (CkCallback &cb)
{
  this->cb = cb;
  nodes[depth].trace(CkCallback(CkReductionTarget(Matrix, doneTrace), thisProxy));
}

/** The reduction target for the trace operation.
 *
 * @param trace The reduction result.
 */
void Matrix::doneTrace (double trace)
{
  this->trace = trace;
  cb.send();
}

/** Get the trace of a matrix.
 *
 * @return The trace.
 */
DoubleMsg * Matrix::getTrace (void)
{
  return new DoubleMsg(trace);
}

/** Add another matrix to this one.
 *
 * @f[ A \leftarrow \alpha A + \beta B @f]
 *
 * where @f$ A @f$ is this Matrix.
 *
 * @param alpha The factor @f$ \alpha @f$.
 * @param beta The factor @f$ \beta @f$.
 * @param B The matrix B.
 * @param cb The callback to call when done.
 */
void Matrix::add (double alpha, double beta, CProxy_Matrix B, CkCallback &cb)
{
  MatrixNodeMsg *BNodes = B.getNodes();
  nodes[depth].add(alpha, beta, BNodes->nodes[depth], CkCallbackResumeThread());
  thisProxy.setNorm(CkCallbackResumeThread());
  cb.send();
}
/** Set this Matrix equal to another matrix.
 *
 * @f[ A \leftarrow B @f]
 *
 * where @f$ A @f$ is this Matrix.
 *
 * @param B The other matrix.
 * @param cb The callback to signal when done.
 */
void Matrix::setEqual (CProxy_Matrix B, CkCallback &cb)
{
  thisProxy.add(0.0, 1.0, B, CkCallbackResumeThread());
  thisProxy.setNorm(CkCallbackResumeThread());
  cb.send();
}

/** Multiply the matrix by a factor.
 *
 * @f[ A \leftarrow \alpha A @f].
 *
 * @param alpha The scale factor.
 * @param cb The callback to signal when done.
 */
void Matrix::scale (double alpha, CkCallback &cb)
{
  DEBUG("scale by %e, depth = %d\n", alpha, depth);
  nodes[depth].scale(alpha, CkCallbackResumeThread());
  thisProxy.setNorm(CkCallbackResumeThread());
  cb.send();
}

/** Add a the scaled identity matrix, i.e. the identity matrix times a scalar
 * to this Matrix.
 *
 * @f[ A \leftarrow \alpha A + \beta I @f]
 *
 * @param alpha The scalar alpha.
 * @param beta The scalar beta.
 * @param cb The callback to signal when done.
 */
void Matrix::addIdentity (double alpha, double beta, CkCallback &cb)
{
  nodes[depth].scale(alpha, CkCallbackResumeThread());
  nodes[depth].addIdentity(beta, CkCallbackResumeThread());
  thisProxy.setNorm(CkCallbackResumeThread());
  cb.send();
}

#include "matrix.def.h"
