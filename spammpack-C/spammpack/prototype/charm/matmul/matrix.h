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

    /** The tree depth of the matrix. */
    int depth;

    /** The padded size of the matrix. */
    int NPadded;

    /** The Nodes per tier. The variable Matrix::nodes points to an array of
     * size Matrix::depth containing the Node proxy on each tier. */
    CProxy_Node *nodes;

    /** The PE map of the nodes at the leaf node tier. */
    int *PEMap;

    /** The norm of the elements in the PEMap. */
    double *PEMap_norm;

    /** A callback. */
    CkCallback cb;

  public:

    Matrix (int N, int blocksize, int nameLength, char *name);
    ~Matrix (void);
    MatrixInfoMsg * info (void);
    DenseMatrixMsg * toDense (void);
    MatrixNodeMsg * getNodes (int tier);
    void updatePEMap (CkCallback &cb);
    void donePEMap (CkReductionMsg *data);
    void set (int N, double *A, CkCallback &cb);
    void setNorm (CkCallback &cb);
};

#endif
