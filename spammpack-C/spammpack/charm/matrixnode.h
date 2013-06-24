#ifndef __MATRIXNODE_
#define __MATRIXNODE_

#include "matrixnode.decl.h"

class MatrixNode : public CBase_MatrixNode
{
  private:

    int NLower[2];
    int NUpper[2];
    int chunksize;
    float *block;

    CProxy_MatrixNode *child[4];

  public:

    MatrixNode (int NLower[2], int NUpper[2], int chunksize);
    NodeInfoMsg * getInfo ();
    BlockMsg * getBlock ();
    void get (int i, int j, CkFuture f);
    void set (int i, int j, float aij, CkFuture f);
    void multiply (CProxy_MatrixNode A, CProxy_MatrixNode B, CkFuture f);
};

#endif
