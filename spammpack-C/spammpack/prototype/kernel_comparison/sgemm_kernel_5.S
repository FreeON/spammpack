#define loops     %rdi
#define N         %rsi
#define alpha     %xmm0
#define A_dilated %rdx
#define B         %rcx
#define C         %r8

#define MULTIPLY

  .text
  .align 256
  .global sgemm_kernel_5
  .type sgemm_kernel_5, @function

  /* Function declaration in C.
   *
   * void
   * sgemm_kernel_5 (const unsigned int loops,
   *     const unsigned int N,
   *     float alpha,
   *     float *restrict A_dilated,
   *     float *restrict B,
   *     float *restrict C);
   */

sgemm_kernel_5:

  /* First we copy the first double word to the other double words in xmm0. */
  pshufd $0x00, alpha, alpha

  /* Loop. */
  movq loops, %rax

loop:

#ifdef MULTIPLY
  movaps 0*4*4(B), %xmm1 /* B_{11} B_{12} B_{13} B_{14} */
  movaps 0*4*4(A_dilated), %xmm12 /* A_{11} A_{11} A_{11} A_{11} */
  mulps %xmm1, %xmm12

  movaps 4*4*4(A_dilated), %xmm13 /* A_{21} A_{21} A_{21} A_{21} */
  mulps %xmm1, %xmm13

  movaps 8*4*4(A_dilated), %xmm14 /* A_{31} A_{31} A_{31} A_{31} */
  mulps %xmm1, %xmm14

  movaps 12*4*4(A_dilated), %xmm15 /* A_{41} A_{41} A_{41} A_{41} */
  mulps %xmm1, %xmm15

  movaps 1*4*4(B), %xmm2 /* B_{21} B_{22} B_{23} B_{24} */
  movaps 1*4*4(A_dilated), %xmm3 /* A_{12} A_{12} A_{12} A_{12} */
  mulps %xmm2, %xmm3

  movaps 5*4*4(A_dilated), %xmm4 /* A_{22} A_{22} A_{22} A_{22} */
  mulps %xmm2, %xmm4

  movaps 9*4*4(A_dilated), %xmm5 /* A_{32} A_{32} A_{32} A_{32} */
  mulps %xmm2, %xmm5

  movaps 13*4*4(A_dilated), %xmm6 /* A_{42} A_{42} A_{42} A_{42} */
  mulps %xmm2, %xmm6

  movaps 2*4*4(B), %xmm7 /* B_{31} B_{32} B_{33} B_{34} */
  movaps 2*4*4(A_dilated), %xmm8 /* A_{13} A_{13} A_{13} A_{13} */
  mulps %xmm7, %xmm8
  addps %xmm3, %xmm12

  movaps 6*4*4(A_dilated), %xmm9 /* A_{23} A_{23} A_{23} A_{23} */
  mulps %xmm7, %xmm9
  addps %xmm4, %xmm13

  movaps 10*4*4(A_dilated), %xmm10 /* A_{33} A_{33} A_{33} A_{33} */
  mulps %xmm7, %xmm10
  addps %xmm5, %xmm14

  movaps 14*4*4(A_dilated), %xmm11 /* A_{43} A_{43} A_{43} A_{43} */
  mulps %xmm7, %xmm11
  addps %xmm6, %xmm15

  movaps 3*4*4(B), %xmm1 /* B_{41} B_{42} B_{43} B_{44} */
  movaps 3*4*4(A_dilated), %xmm2 /* A_{14} A_{14} A_{14} A_{14} */
  mulps %xmm1, %xmm2
  addps %xmm8, %xmm12

  movaps 7*4*4(A_dilated), %xmm3 /* A_{24} A_{24} A_{24} A_{24} */
  mulps %xmm1, %xmm3
  addps %xmm9, %xmm13

  movaps 11*4*4(A_dilated), %xmm4 /* A_{34} A_{34} A_{34} A_{34} */
  mulps %xmm1, %xmm4
  addps %xmm10, %xmm14

  movaps 15*4*4(A_dilated), %xmm5 /* A_{44} A_{44} A_{44} A_{44} */
  mulps %xmm1, %xmm5
  addps %xmm11, %xmm15

  movaps 0*4*4(C), %xmm1
  addps %xmm2, %xmm12
  addps %xmm3, %xmm13

  movaps 1*4*4(C), %xmm2
  addps %xmm4, %xmm14
  addps %xmm5, %xmm15

  movaps 2*4*4(C), %xmm3
  mulps %xmm0, %xmm12
  mulps %xmm0, %xmm13

  movaps 3*4*4(C), %xmm4
  mulps %xmm0, %xmm14
  mulps %xmm0, %xmm15

  addps %xmm1, %xmm12
  addps %xmm2, %xmm13
  movaps %xmm12, 0*4*4(C)

  addps %xmm3, %xmm14
  addps %xmm4, %xmm15
  movaps %xmm13, 1*4*4(C)

  movaps %xmm14, 2*4*4(C)
  movaps %xmm15, 3*4*4(C)
#endif

  sub $1, %rax
  jg loop

  /* We are done. */
  ret

  .size sgemm_kernel_5, .-sgemm_kernel_5
