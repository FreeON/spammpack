#ifndef __MATRIX_H
#define __MATRIX_H

#include "Matrix.decl.h"
#include "Messages.h"

#include <string>

class Matrix : public CBase_Matrix
{
  private:

    int N;
    int blocksize;
    int depth;
    int NPadded;

    CProxy_Node *root;

  public:

    Matrix (int N, int blocksize);
    MatrixMsg * info ();
    DoubleMsg * get (int i, int j);
    EmptyMsg  * set (int i, int j, double aij);
    EmptyMsg  * setBlock (int iLower, int jLower, int iUpper, int jUpper, const double *ABlock);
    IntMsg    * matmul (CProxy_Matrix A, CProxy_Matrix B);
};

#endif
