#ifndef __STRASSENOMP_H
#define __STRASSENOMP_H

class Node
{
  private:

    int blocksize;

    int iLower;
    int iUpper;
    int jLower;
    int jUpper;

    Node *child[4];
    double *data;

  public:

    Node (int blocksize, int iLower, int jLower, int iUpper, int jUpper);
    void set (int i, int j, double aij);
};

class Matrix
{
  private:

    int N;
    int blocksize;
    int depth;
    int NPadded;
    Node *root;

  public:

    Matrix (int N, int blocksize);
    void random ();
    void set (int i, int j, double aij);
};

#endif
