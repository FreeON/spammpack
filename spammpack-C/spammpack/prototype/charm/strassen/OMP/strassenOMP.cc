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

  int c;
  const char *short_options = "hN:b:";
  const option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "N", required_argument, NULL, 'N' },
    { "block", required_argument, NULL, 'b' },
    { NULL, 0, NULL, 0 }
  };

  while((c = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
  {
    switch(c)
    {
      case 'h':
        printf("Usage:\n");
        exit(0);
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'b':
        blocksize = strtol(optarg, NULL, 10);
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

  Matrix A = Matrix(N, blocksize);
  Matrix B = Matrix(N, blocksize);
  Matrix C = Matrix(N, blocksize);

  A.random();
}
