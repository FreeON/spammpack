#include "spamm.h"

#include <stdio.h>
#include <stdlib.h>

extern void sgemm_ (char *, char *, int *, int *, int *, float *, float *, int *, float *, int *, float *, float *, int *);

void
spamm_stream_kernel_C (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream)
{
  unsigned int stream_index;

  struct spamm_data_t *A_data;
  struct spamm_data_t *B_data;
  struct spamm_data_t *C_data;

  float *A_dense;
  float *B_dense;
  float *C_dense;

  int N;
  float beta;

  /* Set some variables. */
  N = SPAMM_N_KERNEL;
  beta = 1.0;

  /* Loop through the stream. */
  for (stream_index = 0; stream_index < number_stream_elements; stream_index++)
  {
    A_data = multiply_stream[stream_index].A;
    B_data = multiply_stream[stream_index].B;
    C_data = multiply_stream[stream_index].C;

    A_dense = A_data->block_dense;
    B_dense = B_data->block_dense;
    C_dense = C_data->block_dense;

    /* Multiply the blocks with an external sgemm() call. We are assuming a
     * Fortran interface, hence sgemm_().
     */
    sgemm_("N", "N", &N, &N, &N, &alpha, A_dense, &N, B_dense, &N, &beta, C_dense, &N);
  }
}
