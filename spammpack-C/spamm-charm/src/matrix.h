/** @file
 *
 * The header file for the Matrix class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __MATRIX_H
#define __MATRIX_H

#include "matrix.decl.h"
#include "node.h"
#include "types.h"

/** A matrix.
 */
class Matrix : public CBase_Matrix
{
  private:

    /** The matrix name. This is currently used to identify the matrix. */
    char *name;

    /** The matrix size. */
    int N;

    /** The submatrix size at the lowest tier. */
    int blocksize;

    /** The size of the basic submatrices. */
    int N_basic;

    /** The tree depth of the matrix. */
    int depth;

    /** The padded size of the matrix. */
    int NPadded;

    /** The trace of the matrix. Needs to be updated with Matrix::updateTrace(). */
    double trace;

    /** The Nodes per tier. The variable Matrix::nodes points to an array of
     * size Matrix::depth containing the Node proxy on each tier. */
    CProxy_Node *nodes;

    /** The PEMap, updated by calling updatePEMap(). */
    int *PEMap;

    /** The norms of the PEMap. */
    double *PEMap_norm;

    /** A callback. */
    CkCallback cb;

  public:

    Matrix (int initialPE, bool alignPEs, int N, int blocksize, int N_basic,
        int nameLength, char *name);
    ~Matrix (void);
    void init (CkCallback &cb);
    MatrixInfoMsg * info (void);
    DenseMatrixMsg * toDense (void);
    MatrixNodeMsg * getNodes (void);
    void updatePEMap (CkCallback &cb);
    void donePEMap (CkReductionMsg *data);
    void set (int N, double *A, CkCallback &cb);
    void setNorm (CkCallback &cb);
    PEMapMsg * getPEMap (void);
    void updateTrace (CkCallback &cb);
    DoubleMsg * getTrace (void);
    void doneTrace (double trace);
    void add (double alpha, double beta, CProxy_Matrix B, CkCallback &cb);
    void setEqual (CProxy_Matrix B, CkCallback &cb);
    void scale (double alpha, CkCallback &cb);
    void addIdentity (double alpha, double beta, CkCallback &cb);
};

#endif
