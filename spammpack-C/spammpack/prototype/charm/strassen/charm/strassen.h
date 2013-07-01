#ifndef __STRASSEN_H
#define __STRASSEN_H

#include "strassen.decl.h"

#include <time.h>
#include <string>

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

class Main : public CBase_Main
{
  public:

    Main (CkArgMsg *msg);
};

class Timer
{
  private:

    std::string name;
    bool isRunning;
    struct timespec startTime;
    struct timespec endTime;

  public:

    Timer (std::string name);
    void start ();
    void stop ();
};

#endif
