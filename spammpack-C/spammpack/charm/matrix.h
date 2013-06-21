#ifndef __MATRIX_H
#define __MATRIX_H

#include "matrix.decl.h"

class Matrix : public CBase_Matrix
{
  private:

    int N;
    int chunksize;
    int NPadded;
    CProxy_MatrixNode *root;

  public:

    Matrix (int N, int chunksize);
    GetMsg* get (int i, int j);
    EmptyMsg* set (int i, int j, float aij);
    EmptyMsg* print ();
};

#endif
