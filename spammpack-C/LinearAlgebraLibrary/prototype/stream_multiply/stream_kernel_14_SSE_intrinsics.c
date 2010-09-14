#include <stdlib.h>
#include <xmmintrin.h>

struct multiply_stream_t
{
  float *A_block;
  float *B_block;
  float *C_block;
  char  mask[8];
};

void
stream_kernel_14_SSE_intrinsics (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream)
{
  short int i;
  unsigned int stream_index;
  unsigned int max_stream_index;

  float *restrict A;
  float *restrict B;
  float *restrict C;

  char *restrict mask;

  __m128 alpha_row;

  __m128 A_element_1;
  __m128 B_row_1;
  __m128 C_row_1;

  __m128 A_element_2;
  __m128 B_row_2;
  __m128 C_row_2;

  __m128 A_element_3;
  __m128 B_row_3;
  __m128 C_row_3;

  __m128 A_element_4;
  __m128 B_row_4;
  __m128 C_row_4;

  __m128 A_element_5;
  __m128 B_row_5;
  __m128 C_row_5;

  __m128 A_element_6;
  __m128 B_row_6;
  __m128 C_row_6;

  __m128 A_element_7;
  __m128 B_row_7;
  __m128 C_row_7;

  __m128 A_element_8;
  __m128 B_row_8;
  __m128 C_row_8;

  /* Divide number of stream elements by 8 to simulate stride of 8. */
  max_stream_index = number_stream_elements/8;

  alpha_row = _mm_set1_ps(alpha);

  for (stream_index = 0; stream_index < max_stream_index; stream_index++)
  {
    /* Load pointers to matrix data blocks. */
    A = multiply_stream[stream_index].A_block;
    B = multiply_stream[stream_index].B_block;
    C = multiply_stream[stream_index].C_block;
    mask = multiply_stream[stream_index].mask;

    if (mask[0]) /* A(1,1)*B(1,1) = C(1,1). */
    {
      for (i = 0; i < 4; i++)
      {
        A_element_1 = _mm_load_ps(&A[(i*4+0)*4]);
        B_row_1 = _mm_load_ps(&B[0*4]);
        C_row_1 = _mm_mul_ps(A_element_1, B_row_1);

        A_element_1 = _mm_load_ps(&A[(i*4+1)*4]);
        B_row_1 = _mm_load_ps(&B[1*4]);
        C_row_1 = _mm_add_ps(_mm_mul_ps(A_element_1, B_row_1), C_row_1);

        A_element_1 = _mm_load_ps(&A[(i*4+2)*4]);
        B_row_1 = _mm_load_ps(&B[2*4]);
        C_row_1 = _mm_add_ps(_mm_mul_ps(A_element_1, B_row_1), C_row_1);

        A_element_1 = _mm_load_ps(&A[(i*4+3)*4]);
        B_row_1 = _mm_load_ps(&B[3*4]);
        C_row_1 = _mm_add_ps(_mm_mul_ps(A_element_1, B_row_1), C_row_1);

        C_row_1 = _mm_mul_ps(alpha_row, C_row_1);
        C_row_1 = _mm_add_ps(_mm_load_ps(&C[i*4]), C_row_1);
        _mm_store_ps(&C[i*4], C_row_1);
      }
    }

    /* Make temporal locality in registers more explicit to the compiler. */

    if (mask[1]) /* A(1,2)*B(2,1) = C(1,1). */
    {
      for (i = 0; i < 4; i++)
      {
        A_element_2 = _mm_load_ps(&A[(i*4+0)*4+64]);
        B_row_2 = _mm_load_ps(&B[0*4+32]);
        C_row_2 = _mm_mul_ps(A_element_2, B_row_2);

        A_element_2 = _mm_load_ps(&A[(i*4+1)*4+64]);
        B_row_2 = _mm_load_ps(&B[1*4+32]);
        C_row_2 = _mm_add_ps(_mm_mul_ps(A_element_2, B_row_2), C_row_2);

        A_element_2 = _mm_load_ps(&A[(i*4+2)*4+64]);
        B_row_2 = _mm_load_ps(&B[2*4+32]);
        C_row_2 = _mm_add_ps(_mm_mul_ps(A_element_2, B_row_2), C_row_2);

        A_element_2 = _mm_load_ps(&A[(i*4+3)*4+64]);
        B_row_2 = _mm_load_ps(&B[3*4+32]);
        C_row_2 = _mm_add_ps(_mm_mul_ps(A_element_2, B_row_2), C_row_2);

        C_row_2 = _mm_mul_ps(alpha_row, C_row_2);
        C_row_2 = _mm_add_ps(_mm_load_ps(&C[i*4]), C_row_2);
        _mm_store_ps(&C[i*4], C_row_2);
      }
    }

    if (mask[2]) /* A(2,1)*B(1,1) = C(2,1). */
    {
      for (i = 0; i < 4; i++)
      {
        A_element_3 = _mm_load_ps(&A[(i*4+0)*4+128]);
        B_row_3 = _mm_load_ps(&B[0*4]);
        C_row_3 = _mm_mul_ps(A_element_3, B_row_3);

        A_element_3 = _mm_load_ps(&A[(i*4+1)*4+128]);
        B_row_3 = _mm_load_ps(&B[1*4]);
        C_row_3 = _mm_add_ps(_mm_mul_ps(A_element_3, B_row_3), C_row_3);

        A_element_3 = _mm_load_ps(&A[(i*4+2)*4+128]);
        B_row_3 = _mm_load_ps(&B[2*4]);
        C_row_3 = _mm_add_ps(_mm_mul_ps(A_element_3, B_row_3), C_row_3);

        A_element_3 = _mm_load_ps(&A[(i*4+3)*4+128]);
        B_row_3 = _mm_load_ps(&B[3*4]);
        C_row_3 = _mm_add_ps(_mm_mul_ps(A_element_3, B_row_3), C_row_3);

        C_row_3 = _mm_mul_ps(alpha_row, C_row_3);
        C_row_3 = _mm_add_ps(_mm_load_ps(&C[i*4+32]), C_row_3);
        _mm_store_ps(&C[i*4+32], C_row_3);
      }
    }

    if (mask[3]) /* A(2,2)*B(2,1) = C(2,1). */
    {
      for (i = 0; i < 4; i++)
      {
        A_element_4 = _mm_load_ps(&A[(i*4+0)*4+192]);
        B_row_4 = _mm_load_ps(&B[0*4+32]);
        C_row_4 = _mm_mul_ps(A_element_4, B_row_4);

        A_element_4 = _mm_load_ps(&A[(i*4+1)*4+192]);
        B_row_4 = _mm_load_ps(&B[1*4+32]);
        C_row_4 = _mm_add_ps(_mm_mul_ps(A_element_4, B_row_4), C_row_4);

        A_element_4 = _mm_load_ps(&A[(i*4+2)*4+192]);
        B_row_4 = _mm_load_ps(&B[2*4+32]);
        C_row_4 = _mm_add_ps(_mm_mul_ps(A_element_4, B_row_4), C_row_4);

        A_element_4 = _mm_load_ps(&A[(i*4+3)*4+192]);
        B_row_4 = _mm_load_ps(&B[3*4+32]);
        C_row_4 = _mm_add_ps(_mm_mul_ps(A_element_4, B_row_4), C_row_4);

        C_row_4 = _mm_mul_ps(alpha_row, C_row_4);
        C_row_4 = _mm_add_ps(_mm_load_ps(&C[i*4+32]), C_row_4);
        _mm_store_ps(&C[i*4+32], C_row_4);
      }
    }

    if (mask[4]) /* A(1,1)*B(1,2) = C(1,2). */
    {
      for (i = 0; i < 4; i++)
      {
        A_element_5 = _mm_load_ps(&A[(i*4+0)*4]);
        B_row_5 = _mm_load_ps(&B[0*4+16]);
        C_row_5 = _mm_mul_ps(A_element_5, B_row_5);

        A_element_5 = _mm_load_ps(&A[(i*4+1)*4]);
        B_row_5 = _mm_load_ps(&B[1*4+16]);
        C_row_5 = _mm_add_ps(_mm_mul_ps(A_element_5, B_row_5), C_row_5);

        A_element_5 = _mm_load_ps(&A[(i*4+2)*4]);
        B_row_5 = _mm_load_ps(&B[2*4+16]);
        C_row_5 = _mm_add_ps(_mm_mul_ps(A_element_5, B_row_5), C_row_5);

        A_element_5 = _mm_load_ps(&A[(i*4+3)*4]);
        B_row_5 = _mm_load_ps(&B[3*4+16]);
        C_row_5 = _mm_add_ps(_mm_mul_ps(A_element_5, B_row_5), C_row_5);

        C_row_5 = _mm_mul_ps(alpha_row, C_row_5);
        C_row_5 = _mm_add_ps(_mm_load_ps(&C[i*4+16]), C_row_5);
        _mm_store_ps(&C[i*4+16], C_row_5);
      }
    }

    if (mask[5]) /* A(1,2)*B(2,2) = C(1,2). */
    {
      for (i = 0; i < 4; i++)
      {
        A_element_6 = _mm_load_ps(&A[(i*4+0)*4+64]);
        B_row_6 = _mm_load_ps(&B[0*4+48]);
        C_row_6 = _mm_mul_ps(A_element_6, B_row_6);

        A_element_6 = _mm_load_ps(&A[(i*4+1)*4+64]);
        B_row_6 = _mm_load_ps(&B[1*4+48]);
        C_row_6 = _mm_add_ps(_mm_mul_ps(A_element_6, B_row_6), C_row_6);

        A_element_6 = _mm_load_ps(&A[(i*4+2)*4+64]);
        B_row_6 = _mm_load_ps(&B[2*4+48]);
        C_row_6 = _mm_add_ps(_mm_mul_ps(A_element_6, B_row_6), C_row_6);

        A_element_6 = _mm_load_ps(&A[(i*4+3)*4+64]);
        B_row_6 = _mm_load_ps(&B[3*4+48]);
        C_row_6 = _mm_add_ps(_mm_mul_ps(A_element_6, B_row_6), C_row_6);

        C_row_6 = _mm_mul_ps(alpha_row, C_row_6);
        C_row_6 = _mm_add_ps(_mm_load_ps(&C[i*4+16]), C_row_6);
        _mm_store_ps(&C[i*4+16], C_row_6);
      }
    }

    if (mask[6]) /* A(2,1)*B(1,2) = C(2,2). */
    {
      for (i = 0; i < 4; i++)
      {
        A_element_7 = _mm_load_ps(&A[(i*4+0)*4+128]);
        B_row_7 = _mm_load_ps(&B[0*4+16]);
        C_row_7 = _mm_mul_ps(A_element_7, B_row_7);

        A_element_7 = _mm_load_ps(&A[(i*4+1)*4+128]);
        B_row_7 = _mm_load_ps(&B[1*4+16]);
        C_row_7 = _mm_add_ps(_mm_mul_ps(A_element_7, B_row_7), C_row_7);

        A_element_7 = _mm_load_ps(&A[(i*4+2)*4+128]);
        B_row_7 = _mm_load_ps(&B[2*4+16]);
        C_row_7 = _mm_add_ps(_mm_mul_ps(A_element_7, B_row_7), C_row_7);

        A_element_7 = _mm_load_ps(&A[(i*4+3)*4+128]);
        B_row_7 = _mm_load_ps(&B[3*4+16]);
        C_row_7 = _mm_add_ps(_mm_mul_ps(A_element_7, B_row_7), C_row_7);

        C_row_7 = _mm_mul_ps(alpha_row, C_row_7);
        C_row_7 = _mm_add_ps(_mm_load_ps(&C[i*4+48]), C_row_7);
        _mm_store_ps(&C[i*4+48], C_row_7);
      }
    }

    if (mask[7]) /* A(2,2)*B(2,2) = C(2,2). */
    {
      for (i = 0; i < 4; i++)
      {
        A_element_8 = _mm_load_ps(&A[(i*4+0)*4+192]);
        B_row_8 = _mm_load_ps(&B[0*4+48]);
        C_row_8 = _mm_mul_ps(A_element_8, B_row_8);

        A_element_8 = _mm_load_ps(&A[(i*4+1)*4+192]);
        B_row_8 = _mm_load_ps(&B[1*4+48]);
        C_row_8 = _mm_add_ps(_mm_mul_ps(A_element_8, B_row_8), C_row_8);

        A_element_8 = _mm_load_ps(&A[(i*4+2)*4+192]);
        B_row_8 = _mm_load_ps(&B[2*4+48]);
        C_row_8 = _mm_add_ps(_mm_mul_ps(A_element_8, B_row_8), C_row_8);

        A_element_8 = _mm_load_ps(&A[(i*4+3)*4+192]);
        B_row_8 = _mm_load_ps(&B[3*4+48]);
        C_row_8 = _mm_add_ps(_mm_mul_ps(A_element_8, B_row_8), C_row_8);

        C_row_8 = _mm_mul_ps(alpha_row, C_row_8);
        C_row_8 = _mm_add_ps(_mm_load_ps(&C[i*4+48]), C_row_8);
        _mm_store_ps(&C[i*4+48], C_row_8);
      }
    }
  }
}
