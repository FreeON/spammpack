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

    MatrixNode (const int NLower[2], const int NUpper[2],
        const int chunksize);
    void get (const int i, const int j, float &aij);
    void set (const int i, const int j, const float aij, CkFuture f);
};

#endif
