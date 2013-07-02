#include "strassen.h"
#include "Messages.h"
#include "Timer.h"
#include "Utilities.h"

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>

#define TOLERANCE 1e-10

/** Multiply two matrices, @f$ C \leftarrow A \times B @f$.
 *
 * @param N The matrix dimension.
 * @Param A Matrix A.
 * @Param B Matrix B.
 *
 * @return The result, matrix C.
 */
double * multiply (int N, double *A, double *B)
{
  double *C = new double[N*N];

  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
    {
      C[i*N+j] = 0;

      for(int k = 0; k < N; k++)
      {
        C[i*N+j] += A[i*N+k]*B[k*N+j];
      }
    }
  }

  return C;
}

/** Convert a matrix in dense row-major storage to a Matrix.
 *
 * @param N The matrix dimension.
 * @param A Pointer to the dense matrix.
 */
void convert (int N, CProxy_Matrix A, double *ADense)
{
  LOG_DEBUG("converting matrix\n");
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
    {
      EmptyMsg *msg = A.set(i, j, ADense[i*N+j]);
      delete msg;
    }
  }
}

/** Set a matrix to zero.
 *
 * @param A The matrix.
 */
void zero (CProxy_Matrix A)
{
  LOG_DEBUG("zeroing matrix\n");
  MatrixMsg *AInfo = A.info();
  for(int i = 0; i < AInfo->N; i++) {
    for(int j = 0; j < AInfo->N; j++)
    {
      EmptyMsg *msg = A.set(i, j, 0.0);
      delete msg;
    }
  }
  delete AInfo;
}

/** Print a matrix.
 *
 * @param name The name of the matrix.
 * @param A The matrix.
 */
void print (std::string name, CProxy_Matrix A)
{
  MatrixMsg *AInfo = A.info();
  std::ostringstream o;
  o.setf(std::ios::fixed);
  o << name << std::endl;
  for(int i = 0; i < AInfo->N; i++) {
    for(int j = 0; j < AInfo->N; j++)
    {
      DoubleMsg *msg = A.get(i, j);
      o << " " << msg->x;
      delete msg;
    }
    o << std::endl;
  }
  CkPrintf("%s", o.str().c_str());
  delete AInfo;
}

/** Allocate and fill a random NxN matrix.
 *
 * @param N The size of the matrix.
 *
 * @return The newly allocated matrix.
 */
double * randomDense (int N)
{
  double *A = new double[N*N];

  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
    {
      A[i*N+j] = rand()/(double) RAND_MAX;
    }
  }

  return A;
}

Main::Main (CkArgMsg *msg)
{
  /* Matrix size. */
  int N = 1;

  /* The submatrix size at the leaf nodes. */
  int blocksize = 1;

  /* Whether to print some debugging stuff. */
  bool debug = false;

  /* Whether to verify the correctness of the matrix product. */
  bool verify = false;

  int c;
  const char *short_options = "hN:b:dv";
  const option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "N", required_argument, NULL, 'N' },
    { "block", required_argument, NULL, 'b' },
    { "debug", no_argument, NULL, 'd' },
    { "verify", no_argument, NULL, 'v' },
    { NULL, 0, NULL, 0 }
  };

  while((c = getopt_long(msg->argc, msg->argv, short_options, long_options, NULL)) != -1)
  {
    switch(c)
    {
      case 'h':
        CkPrintf("Usage: strassen [options]\n");
        CkPrintf("\n");
        CkPrintf("{ --help | -h }       This help\n");
        CkPrintf("-N N                  The matrix size, NxN\n");
        CkPrintf("{ --block | -b } N    The block size, NxN\n");
        CkPrintf("{ --debug | -d }      Print the matrices for debugging\n");
        CkPrintf("{ --verify | -v }     Verify correctness of matrix product\n");
        CkExit();
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'b':
        blocksize = strtol(optarg, NULL, 10);
        break;

      case 'd':
        debug = !debug;
        break;

      case 'v':
        verify = !verify;
        break;

      default:
        CkPrintf("illegal command line argument\n");
        CkExit();
        break;
    }
  }

  if(optind < msg->argc)
  {
    CkPrintf("additional command line arguments, that will be ignored\n");
    for(int i = optind; i < msg->argc; i++)
    {
      CkPrintf("ignoring argument \"%s\"\n", msg->argv[i]);
    }
  }

  if(N < 1 || blocksize < 1)
  {
    CkPrintf("matrix dimensions should be > 0\n");
    CkExit();
  }

  /* Enter threaded main method. */
  thisProxy.run(N, blocksize, debug, verify);
}

void Main::run (int N, int blocksize, bool debug, bool verify)
{
  CProxy_Matrix A = CProxy_Matrix::ckNew(N, blocksize);
  CProxy_Matrix B = CProxy_Matrix::ckNew(N, blocksize);
  CProxy_Matrix C = CProxy_Matrix::ckNew(N, blocksize);

  LOG_DEBUG("N = %d, blocksize = %d\n", N, blocksize);
  double *ADense = randomDense(N);
  double *BDense = randomDense(N);

  LOG_DEBUG("converting dense matrices to Matrix\n");
  convert(N, A, ADense);
  convert(N, B, BDense);

  LOG_DEBUG("setting C to zero\n");
  zero(C);

  if(debug)
  {
    print("A:", A);
    print("B:", B);
  }

  Timer timer("multiply");
  timer.start();
  EmptyMsg *msg = C.matmul(A, B); delete msg;
  timer.stop();

  if(debug) print("C:", C);

  if(verify)
  {
    double *CDense = multiply(N, ADense, BDense);
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++)
      {
        DoubleMsg *cij = C.get(i, j);
        if(fabs(CDense[i*N+j]-cij->x) > TOLERANCE)
        {
          CkPrintf("comparison failed for C(%d,%d): %e <-> %e\n",
              i, j, CDense[i*N+j], cij->x);
          CkExit();
        }
        delete cij;
      }
    }
    CkPrintf("result verified\n");
  }

  CkExit();
}

#include "strassen.def.h"
