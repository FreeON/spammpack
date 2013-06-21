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

    Matrix ();
    Matrix (int N, int chunksize);
    EmptyMsg * initialize (int N, int chunksize);
    EmptyMsg * remove();
    IntMsg * getN ();
    IntMsg * getChunksize ();
    FloatMsg * get (int i, int j);
    EmptyMsg * set (int i, int j, float aij);
    EmptyMsg * print ();
};

#endif
