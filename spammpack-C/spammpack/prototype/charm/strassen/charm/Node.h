#ifndef __NODE_H
#define __NODE_H

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

    int blockIndex (int i, int j);

  public:

    Node (int blocksize, int iLower, int jLower, int iUpper, int jUpper);
    void set (int i, int j, double aij);
    double get (int i, int j);
    void matmul (Node A, Node B);
};

#endif
