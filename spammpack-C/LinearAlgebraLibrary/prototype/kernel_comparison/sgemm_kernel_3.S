#define N         %rdi
#define alpha     %xmm0
#define A_dilated %rsi
#define B         %rdx
#define C         %rcx

  .text
  .align 256
  .global sgemm_kernel_3
  .type sgemm_kernel_3, @function

  /* Function declaration in C.
   *
   * void
   * sgemm_kernel_3 (const unsigned int N,
   *     float alpha,
   *     float *restrict A_dilated,
   *     float *restrict B,
   *     float *restrict C);
   */

sgemm_kernel_3:

  /* Create some space on the stack. */
  movq %rsp, %rax

  /* First we copy the first double word to the other double words in xmm0. */
  pshufd $0x00, alpha, alpha

  /* Zero out result matrix. */
  pxor %xmm8,  %xmm8
  pxor %xmm9,  %xmm9
  pxor %xmm10, %xmm10
  pxor %xmm11, %xmm11

  /* Load elements. */
  movaps  0*4*4(B),         %xmm1 /* B_{11} B_{12} B_{13} B_{14} */
  movaps  0*4*4(A_dilated), %xmm2 /* A_{11} A_{11} A_{11} A_{11} */
  movaps  4*4*4(A_dilated), %xmm3 /* A_{21} A_{21} A_{21} A_{21} */
  movaps  8*4*4(A_dilated), %xmm4 /* A_{31} A_{31} A_{31} A_{31} */
  movaps 12*4*4(A_dilated), %xmm5 /* A_{41} A_{41} A_{41} A_{41} */

  /* Multiply A*B and store in C. */
  mulps %xmm1, %xmm2
  addps %xmm2, %xmm8

  mulps %xmm1, %xmm3
  addps %xmm3, %xmm9

  mulps %xmm1, %xmm4
  addps %xmm4, %xmm10

  mulps %xmm1, %xmm5
  addps %xmm5, %xmm11

  /* Load elements. */
  movaps  1*4*4(B),         %xmm1 /* B_{21} B_{22} B_{23} B_{24} */
  movaps  1*4*4(A_dilated), %xmm2 /* A_{12} A_{12} A_{12} A_{12} */
  movaps  5*4*4(A_dilated), %xmm3 /* A_{22} A_{22} A_{22} A_{22} */
  movaps  9*4*4(A_dilated), %xmm4 /* A_{32} A_{32} A_{32} A_{32} */
  movaps 13*4*4(A_dilated), %xmm5 /* A_{42} A_{42} A_{42} A_{42} */

  /* Multiply A*B and store in C. */
  mulps %xmm1, %xmm2
  addps %xmm2, %xmm8

  mulps %xmm1, %xmm3
  addps %xmm3, %xmm9

  mulps %xmm1, %xmm4
  addps %xmm4, %xmm10

  mulps %xmm1, %xmm5
  addps %xmm5, %xmm11

  /* Load elements. */
  movaps  2*4*4(B),         %xmm1 /* B_{31} B_{32} B_{33} B_{34} */
  movaps  2*4*4(A_dilated), %xmm2 /* A_{13} A_{13} A_{13} A_{13} */
  movaps  6*4*4(A_dilated), %xmm3 /* A_{23} A_{23} A_{23} A_{23} */
  movaps 10*4*4(A_dilated), %xmm4 /* A_{33} A_{33} A_{33} A_{33} */
  movaps 14*4*4(A_dilated), %xmm5 /* A_{43} A_{43} A_{43} A_{43} */

  /* Multiply A*B and store in C. */
  mulps %xmm1, %xmm2
  addps %xmm2, %xmm8

  mulps %xmm1, %xmm3
  addps %xmm3, %xmm9

  mulps %xmm1, %xmm4
  addps %xmm4, %xmm10

  mulps %xmm1, %xmm5
  addps %xmm5, %xmm11

  /* Load elements. */
  movaps  3*4*4(B),         %xmm1 /* B_{41} B_{42} B_{43} B_{44} */
  movaps  3*4*4(A_dilated), %xmm2 /* A_{14} A_{14} A_{14} A_{14} */
  movaps  7*4*4(A_dilated), %xmm3 /* A_{24} A_{24} A_{24} A_{24} */
  movaps 11*4*4(A_dilated), %xmm4 /* A_{34} A_{34} A_{34} A_{34} */
  movaps 15*4*4(A_dilated), %xmm5 /* A_{44} A_{44} A_{44} A_{44} */

  /* Multiply A*B and store in C. */
  mulps %xmm1, %xmm2
  addps %xmm2, %xmm8

  mulps %xmm1, %xmm3
  addps %xmm3, %xmm9

  mulps %xmm1, %xmm4
  addps %xmm4, %xmm10

  mulps %xmm1, %xmm5
  addps %xmm5, %xmm11

  /* Multiply alpha. */
  mulps %xmm0, %xmm8
  mulps %xmm0, %xmm9
  mulps %xmm0, %xmm10
  mulps %xmm0, %xmm11

  /* Add C from registers back into memory. */
  movaps 0*4*4(C), %xmm0
  addps %xmm0, %xmm8

  movaps 1*4*4(C), %xmm1
  addps %xmm1, %xmm9

  movaps 2*4*4(C), %xmm2
  addps %xmm2, %xmm10

  movaps 3*4*4(C), %xmm3
  addps %xmm3, %xmm11

  movaps %xmm8,  0*4*4(C)
  movaps %xmm9,  1*4*4(C)
  movaps %xmm10, 2*4*4(C)
  movaps %xmm11, 3*4*4(C)

  /* We are done. */
  ret

  .size sgemm_kernel_3, .-sgemm_kernel_3
