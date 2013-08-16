/* @file
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

/** A matrix.
 */
class Matrix : public CBase_Matrix
{
  private:

    /** The matrix size. */
    int N;

    /** The submatrix size at the lowest tier. */
    int blocksize;

    /** The tree depth of the matrix. */
    int depth;

    /** The padded size of the matrix. */
    int NPadded;

    /** The Node matrix at tier == depth. */
    CProxy_Node nodes;

  public:

    Matrix (int N, int blocksize);
    MatrixInfoMsg * info (void);
    DenseMatrixMsg * toDense (void);
};

#endif
