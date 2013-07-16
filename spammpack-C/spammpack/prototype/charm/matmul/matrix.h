#ifndef __MATRIX_H
#define __MATRIX_H

#include "matrix.decl.h"

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
};

#endif
