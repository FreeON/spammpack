#include <stdio.h>
#include <stdlib.h>

#define RANDOM_MATRIX
//#define PRINT_DEBUG

#if defined(SSE)
void
matrix_multiply_SSE (const unsigned int N, float *A, float *B, float *C);
#elif defined(SSE4_1)
void
matrix_multiply_SSE4_1 (const unsigned int N, float *A, float *B, float *C);
#endif

int
main (int argc, char **argv)
{
  float __attribute__ ((aligned (64))) A[4][4];
  float __attribute__ ((aligned (64))) A_dilated[4][4][4];
  float __attribute__ ((aligned (64))) B[4][4];
  float __attribute__ ((aligned (64))) B_transpose[4][4];
  float __attribute__ ((aligned (64))) C[4][4];

  short i, j;

  unsigned int max_N = 1;

  /* Parse command line. */
  if (argc == 2)
  {
    max_N = strtol(argv[1], NULL, 10);
  }

  /* Fill matrix with some random stuff. */
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++)
    {
#ifndef RANDOM_MATRIX
      A[i][j] = i*4+j;
      B[i][j] = i*4+j;
      C[i][j] = i*4+j;
#else
      A[i][j] = rand()/(float) RAND_MAX;
      B[i][j] = rand()/(float) RAND_MAX;
      C[i][j] = rand()/(float) RAND_MAX;
#endif
      B_transpose[j][i] = B[i][j];
      A_dilated[i][j][0] = A[i][j];
      A_dilated[i][j][1] = A[i][j];
      A_dilated[i][j][2] = A[i][j];
      A_dilated[i][j][3] = A[i][j];
    }
  }

#ifdef SSE
  matrix_multiply_SSE(max_N, (float*) &A_dilated[0][0], (float*) &B[0][0], (float*) &C[0][0]);
#elif defined(SSE4_1)
  matrix_multiply_SSE4_1(max_N, (float*) &A[0][0], (float*) &B_transpose[0][0], (float*) &C[0][0]);
#endif

#ifdef PRINT_DEBUG
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++)
    {
      //printf(" %i", (int) C[i][j]);
      printf(" %f", C[i][j]);
    }
    printf("\n");
  }
#endif

  return 0;
}
