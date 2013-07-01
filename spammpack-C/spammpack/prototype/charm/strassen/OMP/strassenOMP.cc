#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include "strassenOMP.h"

int main (int argc, char **argv)
{
  /* Matrix size. */
  int N = 1;

  /* The submatrix size at the leaf nodes. */
  int blocksize = 1;

  /* Whether to print some debugging stuff. */
  bool debug = false;

  /* Whether to verify the correctness of the matrix product. */
  bool verify = false;

  double *ADense;
  double *BDense;

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

  while((c = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
  {
    switch(c)
    {
      case 'h':
        printf("Usage: strassenOMP [options]\n");
        printf("\n");
        printf("{ --help | -h }       This help\n");
        printf("-N N                  The matrix size, NxN\n");
        printf("{ --block | -b } N    The block size, NxN\n");
        printf("{ --debug | -d }      Print the matrices for debugging\n");
        printf("{ --verify | -v }     Verify correctness of matrix product\n");
        exit(0);
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
        printf("illegal command line argument\n");
        exit(1);
        break;
    }
  }

  if(optind < argc)
  {
    printf("additional command line arguments, that will be ignored\n");
  }

  if(N < 1 || blocksize < 1)
  {
    printf("matrix dimensions should be > 0\n");
    exit(1);
  }

#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp master
    printf("running on % 2d cores... ", omp_get_num_threads()); fflush(stdout);
  }
#endif

  Matrix A = Matrix(N, blocksize);
  Matrix B = Matrix(N, blocksize);
  Matrix C = Matrix(N, blocksize);

  A.random();
  B.random();
  C.zero();

  if(debug)
  {
    A.print("A:");
    B.print("B:");
  }

  Timer timer("multiply");
  timer.start();
  C.matmul(A, B);
  timer.stop();

  if(debug) C.print("C:");
}
