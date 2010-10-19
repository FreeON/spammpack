#include <stdlib.h>
#include <xmmintrin.h>

struct multiply_stream_t
{
  float *A_block;
  float *B_block;
  float *C_block;
  float  norm[32];
};

void
stream_kernel_19_SSE_intrinsics (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct multiply_stream_t *multiply_stream)
{

  static const int A_OFFSET[4][4] =
  {
    {
      0, /* (0*4+0)*64 */
      64, /* (0*4+1)*64 */
      128, /* (0*4+2)*64 */
      192 /* (0*4+3)*64 */
    },
    {
      256, /* (1*4+0)*64 */
      320, /* (1*4+1)*64 */
      384, /* (1*4+2)*64 */
      448 /* (1*4+3)*64 */
    },
    {
      512, /* (2*4+0)*64 */
      576, /* (2*4+1)*64 */
      640, /* (2*4+2)*64 */
      704 /* (2*4+3)*64 */
    },
    {
      768, /* (3*4+0)*64 */
      832, /* (3*4+1)*64 */
      896, /* (3*4+2)*64 */
      960 /* (3*4+3)*64 */
    }
  };

  static const int B_OFFSET[4][4] =
  {
    {
      0, /* (0*4+0)*16 */
      16, /* (0*4+1)*16 */
      32, /* (0*4+2)*16 */
      48 /* (0*4+3)*16 */
    },
    {
      64, /* (1*4+0)*16 */
      80, /* (1*4+1)*16 */
      96, /* (1*4+2)*16 */
      112 /* (1*4+3)*16 */
    },
    {
      128, /* (2*4+0)*16 */
      144, /* (2*4+1)*16 */
      160, /* (2*4+2)*16 */
      176 /* (2*4+3)*16 */
    },
    {
      192, /* (3*4+0)*16 */
      208, /* (3*4+1)*16 */
      224, /* (3*4+2)*16 */
      240 /* (3*4+3)*16 */
    }
  };

  static const int C_OFFSET[4][4] =
  {
    {
      0, /* (0*4+0)*16 */
      16, /* (0*4+1)*16 */
      32, /* (0*4+2)*16 */
      48 /* (0*4+3)*16 */
    },
    {
      64, /* (1*4+0)*16 */
      80, /* (1*4+1)*16 */
      96, /* (1*4+2)*16 */
      112 /* (1*4+3)*16 */
    },
    {
      128, /* (2*4+0)*16 */
      144, /* (2*4+1)*16 */
      160, /* (2*4+2)*16 */
      176 /* (2*4+3)*16 */
    },
    {
      192, /* (3*4+0)*16 */
      208, /* (3*4+1)*16 */
      224, /* (3*4+2)*16 */
      240 /* (3*4+3)*16 */
    }
  };

  short int i, j, k, l;
  unsigned int stream_index;
  unsigned int max_stream_index;

  float *restrict A;
  float *restrict B;
  float *restrict C;

  float *restrict norm;

  __m128 alpha_row;

  __m128 A_element;
  __m128 B_row;
  __m128 C_row[4];

  /* Divide number of stream elements by 64 to simulate stride of 64. */
  max_stream_index = number_stream_elements/64;

  alpha_row = _mm_set1_ps(alpha);

  for (stream_index = 0; stream_index < max_stream_index; stream_index++)
  {
    /* Load pointers to matrix data blocks. */
    A = multiply_stream[stream_index].A_block;
    B = multiply_stream[stream_index].B_block;
    C = multiply_stream[stream_index].C_block;
    norm = multiply_stream[stream_index].norm;

    for (i = 0; i < 4; i++) {
      for (j = 0; j < 4; j++)
      {
        C_row[0] = _mm_setzero_ps();
        C_row[1] = _mm_setzero_ps();
        C_row[2] = _mm_setzero_ps();
        C_row[3] = _mm_setzero_ps();

        for (k = 0; k < 4; k++)
        {
          if (norm[i*4+k]*norm[k*4+j+16] < tolerance) { continue; }
          else
          {
            for (l = 0; l < 4; l++)
            {
              A_element = _mm_load_ps(&A[(l*4+0)*4+A_OFFSET[i][k]]);
              B_row = _mm_load_ps(&B[0*4+B_OFFSET[k][j]]);
              C_row[l] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[l]);

              A_element = _mm_load_ps(&A[(l*4+1)*4+A_OFFSET[i][k]]);
              B_row = _mm_load_ps(&B[1*4+B_OFFSET[k][j]]);
              C_row[l] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[l]);

              A_element = _mm_load_ps(&A[(l*4+2)*4+A_OFFSET[i][k]]);
              B_row = _mm_load_ps(&B[2*4+B_OFFSET[k][j]]);
              C_row[l] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[l]);

              A_element = _mm_load_ps(&A[(l*4+3)*4+A_OFFSET[i][k]]);
              B_row = _mm_load_ps(&B[3*4+B_OFFSET[k][j]]);
              C_row[l] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[l]);
            }
          }
        }

        /* Store C block. */
        for (l = 0; l < 4; l++)
        {
          C_row[l] = _mm_mul_ps(alpha_row, C_row[l]);
          C_row[l] = _mm_add_ps(_mm_load_ps(&C[l*4+C_OFFSET[i][j]]), C_row[l]);
          _mm_store_ps(&C[l*4+C_OFFSET[i][j]], C_row[l]);
        }
      }
    }
  }
}
