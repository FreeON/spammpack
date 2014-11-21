#include <stdlib.h>
#include <stdio.h>

void convolute (unsigned int *A_index_array, float *A_norm,
    unsigned int *B_index_array, float *B_norm);

void convolute_C (unsigned int *A_index_array, float *A_norm,
    unsigned int *B_index_array, float *B_norm);

int
main ()
{
  int N = 10;
  int result;
  unsigned int *A_index_array;
  unsigned int *B_index_array;
  float * A_norm;
  float * B_norm;

  A_index_array = malloc(N*sizeof(unsigned int));
  A_index_array[0] = 1;
  A_index_array[1] = 2;
  A_index_array[2] = 3;
  A_index_array[3] = 4;

  A_norm = malloc(N*sizeof(float));
  A_norm[0] = 0.5;
  A_norm[1] = 0.6;
  A_norm[2] = 0.7;
  A_norm[3] = 0.8;

  result = posix_memalign((void**) &B_index_array, 16, N*sizeof(unsigned int));
  if (result != 0)
  {
    printf("error allocating aligned memory\n");
    exit(1);
  }
  B_index_array[0] = 10;
  B_index_array[1] = 20;
  B_index_array[2] = 30;
  B_index_array[3] = 40;

  result = posix_memalign((void**) &B_norm, 16, N*sizeof(float));
  if (result != 0)
  {
    printf("error allocating aligned memory\n");
    exit(1);
  }
  B_norm[0] = 1.5;
  B_norm[1] = 1.6;
  B_norm[2] = 1.7;
  B_norm[3] = 1.8;
  B_norm[4] = 1.9;
  B_norm[5] = 1.95;

  convolute(A_index_array, A_norm, B_index_array, B_norm);

  convolute_C(A_index_array, A_norm, B_index_array, B_norm);
}
