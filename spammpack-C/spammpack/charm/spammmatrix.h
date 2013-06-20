#ifndef __SPAMMMATRIX_H
#define __SPAMMMATRIX_H

#include "spammmatrix.decl.h"

class SpAMMMatrix : public CBase_SpAMMMatrix
{
  private:

    int N;
    int chunksize;
    int NPadded;
    CProxy_MatrixNode *root;

  public:

    SpAMMMatrix (const int N, const int chunksize);
    void get (const int i, const int j, float &aij);
    void set (const int i, const int j, const float aij);
    void print ();
};

#endif
