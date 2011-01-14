#include "spamm.h"
#include "config.h"

#ifdef HAVE_SSE
#include <xmmintrin.h>
#else
#include <stdio.h>
#include <stdlib.h>
#endif

void
spamm_stream_kernel_C (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream)
{
#ifdef HAVE_SSE
  unsigned int stream_index;

  short int i, j, k, l, m;

  struct spamm_data_t *A_data;
  struct spamm_data_t *B_data;
  struct spamm_data_t *C_data;

  float *A;
  float *B;
  float *C;

  float *norm_A;
  float *norm_B;

  unsigned int A_offset;
  unsigned int B_offset;
  unsigned int C_offset;

  __m128 alpha_row;
  __m128 C_row[4];
  __m128 A_element;
  __m128 B_row;

  alpha_row = _mm_set1_ps(alpha);

  for (stream_index = 0; stream_index < number_stream_elements; stream_index++)
  {
    A_data = multiply_stream[stream_index].A;
    B_data = multiply_stream[stream_index].B;
    C_data = multiply_stream[stream_index].C;

    A = A_data->block_dense_dilated;
    B = B_data->block_dense;
    C = C_data->block_dense;

    norm_A = A_data->norm;
    norm_B = B_data->norm;

    for (i = 0; i < SPAMM_N_KERNEL_BLOCK; i++) {
      for (j = 0; j < SPAMM_N_KERNEL_BLOCK; j++)
      {
        /* Zero C block. */
        for (l = 0; l < 4; l++)
        {
          C_row[l] = _mm_setzero_ps();
        }

        /* Calculate C offset. */
        C_offset = i*SPAMM_N_KERNEL_BLOCK+j;

        for (k = 0; k < SPAMM_N_KERNEL_BLOCK; k += SPAMM_N_STRIDE)
        {
          for (m = 0; m < SPAMM_N_STRIDE; m++)
          {
            /* Calculate offsets. */
            A_offset = i*SPAMM_N_KERNEL_BLOCK+(k+m);
            B_offset = (k+m)*SPAMM_N_KERNEL_BLOCK+j;

            /* Check the norm. */
            if (norm_A[A_offset]*norm_B[B_offset] <= tolerance)
            {
              continue;
            }

            /* Multiply blocks. */
            for (l = 0; l < 4; l++)
            {
              A_element = _mm_load_ps(&A[(l*4+0)*4+A_offset*SPAMM_N_BLOCK*SPAMM_N_BLOCK*4]);
              B_row = _mm_load_ps(&B[0*4+B_offset*SPAMM_N_BLOCK*SPAMM_N_BLOCK]);
              C_row[l] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[l]);

              A_element = _mm_load_ps(&A[(l*4+1)*4+A_offset*SPAMM_N_BLOCK*SPAMM_N_BLOCK*4]);
              B_row = _mm_load_ps(&B[1*4+B_offset*SPAMM_N_BLOCK*SPAMM_N_BLOCK]);
              C_row[l] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[l]);

              A_element = _mm_load_ps(&A[(l*4+2)*4+A_offset*SPAMM_N_BLOCK*SPAMM_N_BLOCK*4]);
              B_row = _mm_load_ps(&B[2*4+B_offset*SPAMM_N_BLOCK*SPAMM_N_BLOCK]);
              C_row[l] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[l]);

              A_element = _mm_load_ps(&A[(l*4+3)*4+A_offset*SPAMM_N_BLOCK*SPAMM_N_BLOCK*4]);
              B_row = _mm_load_ps(&B[3*4+B_offset*SPAMM_N_BLOCK*SPAMM_N_BLOCK]);
              C_row[l] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[l]);
            }
          }
        }

        /* Store C block. */
        for (l = 0; l < 4; l++)
        {
          C_row[l] = _mm_mul_ps(alpha_row, C_row[l]);
          C_row[l] = _mm_add_ps(_mm_load_ps(&C[l*4+C_offset*SPAMM_N_BLOCK*SPAMM_N_BLOCK]), C_row[l]);
          _mm_store_ps(&C[l*4+C_offset*SPAMM_N_BLOCK*SPAMM_N_BLOCK], C_row[l]);
        }
      }
    }
  }
#else
  printf("no SSE, this kernel will not work.\n");
  exit(1);
#endif
}
