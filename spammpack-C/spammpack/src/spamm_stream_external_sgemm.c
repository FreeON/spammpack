#include "spamm.h"

#include <stdio.h>
#include <stdlib.h>

extern void sgemm_ (char *, char *, int *, int *, int *, float *, float *, int *, float *, int *, float *, float *, int *);

void
spamm_stream_external_sgemm (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream,
    const short call_external_sgemm)
{
  unsigned int stream_index;

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
    /* Multiply the blocks with an external sgemm() call. We are assuming a
     * Fortran interface, hence sgemm_().
     */
    if (call_external_sgemm == 1)
    {
      sgemm_("N", "N", &N, &N, &N, &alpha,
          multiply_stream[stream_index].A->block_dense, &N,
          multiply_stream[stream_index].B->block_dense, &N, &beta,
          multiply_stream[stream_index].C->block_dense, &N);
    }
  }
}
