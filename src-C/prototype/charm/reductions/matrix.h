#ifndef __MATRIX_H
#define __MATRIX_H

#include "matrix.decl.h"

class Matrix : public CBase_Matrix
{
  private:

    CkCallback normCB;
    CProxy_Node *root;

  public:

    Matrix (int depth, int childsize);
    void norm (const CkCallback &cb);
    void normDone (double norm);
};

#endif
