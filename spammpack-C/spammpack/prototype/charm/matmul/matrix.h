#ifndef __MATRIX_H
#define __MATRIX_H

#include "matrix.decl.h"
#include "messages.h"
#include "types.h"

class Matrix : public CBase_Matrix
{
  private:

    int N;
    int blocksize;
    int depth;
    int NPadded;

    bool rootNull;
    CProxy_Node root;

  public:

    Matrix (int N, int blocksize);
    MatrixInfoMsg * info ();
    DenseMatrixMsg * getDense ();
    void random (CkCallback &cb);
    void zero (CkCallback &cb);
    void initialize (enum init_t initType, CkCallback &cb);
    void print (CkCallback &cb);
    void multiply (CProxy_Matrix A, CProxy_Matrix B, CkCallback &cb);
};

#endif
