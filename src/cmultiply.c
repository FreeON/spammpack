/** Simple C version of 16x16 single precision matrix multiply.
 *
 * @param C Pointer to dense matrix C.
 * @param A Pointer to dense matrix A.
 * @param B Pointer to dense matrix B.
 */
void
cmultiply_ (float *restrict C, float *restrict A, float *restrict B)
{
  __assume_aligned (A, 64);
  __assume_aligned (B, 64);
  __assume_aligned (C, 64);

  int i, j, k;

  #pragma ivdep
  #pragma vector aligned
  #pragma vector always
  for (i = 0; i < 16; i++) {
    for (j = 0; j < 16; j++) {
      for (k = 0; k < 16; k++)
      {
        C[i+j*16] += A[i+k*16]*B[k+j*16];
      }
    }
  }
}
