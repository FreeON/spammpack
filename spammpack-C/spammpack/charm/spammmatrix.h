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
    void set (const int i, const int j, const float aij);
};

#endif
