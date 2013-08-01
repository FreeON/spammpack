/* @file
 *
 * The header file for the Matrix class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __MATRIX_H
#define __MATRIX_H

#include "messages.h"
#include "types.h"

#include "matrix.decl.h"

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

    /** The array of Node objects. There is one array per tier. */
    CProxy_Node *tierNode;

  public:

    Matrix (int N, int blocksize);
    DenseMatrixMsg * getDense ();
    MatrixInfoMsg * info (int tier);
    void random (CkCallback &cb);
    void zero (CkCallback &cb);
    void decay (double decayConstant, CkCallback &cb);
    void initialize (int initType, double decayConstant);
    void print (CkCallback &cb);
    void printLeafPes (CkCallback &cb);
};

#endif
