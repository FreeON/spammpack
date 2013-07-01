#ifndef __MATRIX_H
#define __MATRIX_H

#include "Matrix.decl.h"

#include <string>

class Matrix : public CBase_Matrix
{
  private:

    int N;
    int blocksize;
    int depth;
    int NPadded;

    Node *root;

  public:

    Matrix (int N, int blocksize);
    void convert (int N, double *A);
    void random ();
    void zero ();
    void set (int i, int j, double aij);
    double get (int i, int j);
    void print (std::string name);
    void matmul (Matrix A, Matrix B);
};

#endif
