#include <xmmintrin.h>

void convolute_C (unsigned int *A_index_array, float *A_norm,
    unsigned int *B_index_array, float *B_norm)
{
  __m128i A_index_xmm;
  __m128 A_norm_xmm;

  A_norm_xmm = _mm_load_ss(&A_norm[0]);
  A_norm_xmm = _mm_shuffle_ps(A_norm_xmm, A_norm_xmm, 0x00);

  A_index_xmm = _mm_cvtsi32_si128(A_index_array[0]);
  A_index_xmm = (__m128i) _mm_shuffle_ps(A_index_xmm, A_index_xmm, 0x00);
}
