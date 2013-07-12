#ifndef __MATRIX_H
#define __MATRIX_H

#include "matrix.decl.h"

class Matrix : public CBase_Matrix
{
  private:

    CProxy_Node *root;

  public:

    Matrix (int depth, int childsize);
    void norm ();
};

#endif
