#include "spamm.h"

/* Compile this file with:
 *
 * gcc -S -I../src kernel.c
 *
 * to get the gcc assembly.
 */

/* The stream elements contain a matrix of NxN basic matrix blocks. */
#ifndef N
#define N 1
#endif

/* The stride to go through a stripe of basic matrix blocks. */
#ifndef N_STRIPE
#define N_STRIPE 1
#endif

#define OFFSET_NORM 24
#define OFFSET_BLOCK_DENSE 192
#define OFFSET_BLOCK_DENSE_DILATED 1216

void
block_multiply (float *C, float *A, float *B)
{
}

void
zero_block (float A[4][4])
{
}

void
norm_check (float norm_A, float norm_B)
{
}

void
spamm_stream_kernel (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream)
{
  short i, j, k;
  unsigned int stream_index;

  struct spamm_data_t *A;
  struct spamm_data_t *B;
  struct spamm_data_t *C;

  float C_acc[4][4];

  for (stream_index = 0; stream_index < number_stream_elements; stream_index++)
  {
    A = multiply_stream[stream_index].A;
    B = multiply_stream[stream_index].B;
    C = multiply_stream[stream_index].C;

    /* Loop over stripe of A and B and accumulate in block of C. */
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        zero_block(C_acc);
        for (k = 0; k < N; k += N_STRIPE)
        {
          /* Check norm. */
          norm_check(*((float*) ((char*) A+OFFSET_NORM)+i),
              *((float*) ((char*) B+OFFSET_NORM)+j));

          /* Multiply blocks. */
          block_multiply((float*) C_acc,
              (float*) ((char*) A+OFFSET_BLOCK_DENSE_DILATED)+i*N*4+k*4,
              (float*) ((char*) B+OFFSET_BLOCK_DENSE)+k*N*4+j*4);
        }
      }
    }
  }
}
