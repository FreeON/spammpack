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
    void random ();
    EmptyMsg * set (int i, int j, double aij);
    DoubleMsg * get (int i, int j);
    void print (std::string name);
    EmptyMsg * matmul (CProxy_Matrix A, CProxy_Matrix B);
    MatrixMsg * info ();
};

#endif
