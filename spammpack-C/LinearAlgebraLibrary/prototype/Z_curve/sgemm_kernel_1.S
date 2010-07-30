#define alpha     %xmm0
#define A         %rdi
#define B         %rsi
#define beta      %xmm1
#define C         %rdx
#define OLD_STACK %r11

/* Function interface:
 *
 * void
 * sgemm_kernel_1 (const float alpha,
 *     const float *restrict A,
 *     const float *restrict B,
 *     const float beta,
 *     float *restrict C);
 */

  .text
  .align 256
  .global sgemm_kernel_1
  .type sgemm_kernel_1, @function

sgemm_kernel_1:

  /* Copy A into buffer. We will copy each matrix element of A 4 times
   * and arrange the elements in the buffer according to the following
   * storage layout:
   *
   * A_{11} A_{12} A_{13} A_{14}
   * A_{21} A_{22} A_{23} ...
   *
   * Copied to buffer as:
   *
   * A_{11} A_{11} A_{11} A_{11}
   * A_{12} A_{12} A_{12} A_{12}
   * A_{13} A_{13} ...
   */

  /* Make space for the buffer for B on the stack. */
  movq %rsp, OLD_STACK /* Save old stack. */
  subq $(64 * 4), %rsp
  andq $-4096, %rsp /* Align stack on a page boundary. */

  /* First we copy the first double word to the other double words in xmm0 and
   * xmm1. */
  pshufd $0x00, alpha, alpha
  //pshufd $0x00, beta, beta

  /* Apply beta to C. */
  //movaps 0*4*4(C), %xmm2
  //movaps 1*4*4(C), %xmm3
  //movaps 2*4*4(C), %xmm4
  //movaps 3*4*4(C), %xmm5

  //mulps beta, %xmm2
  //mulps beta, %xmm3
  //mulps beta, %xmm4
  //mulps beta, %xmm5

  //movaps %xmm2, 0*4*4(C)
  //movaps %xmm3, 1*4*4(C)
  //movaps %xmm4, 2*4*4(C)
  //movaps %xmm5, 3*4*4(C)

  /* Copy matrix elements of A into SSE registers. */
  movaps 0*4*4(A), %xmm3
  movaps 1*4*4(A), %xmm7
  movaps 2*4*4(A), %xmm11
  movaps 3*4*4(A), %xmm15

  /* Apply alpha. */
  mulps alpha, %xmm3
  mulps alpha, %xmm7
  mulps alpha, %xmm11
  mulps alpha, %xmm15

  /* Dilate A. */
  pshufd $0x00, %xmm3, %xmm0
  pshufd $0x55, %xmm3, %xmm1
  pshufd $0xaa, %xmm3, %xmm2
  pshufd $0xff, %xmm3, %xmm3

  pshufd $0x00, %xmm7, %xmm4
  pshufd $0x55, %xmm7, %xmm5
  pshufd $0xaa, %xmm7, %xmm6
  pshufd $0xff, %xmm7, %xmm7

  pshufd $0x00, %xmm11, %xmm8
  pshufd $0x55, %xmm11, %xmm9
  pshufd $0xaa, %xmm11, %xmm10
  pshufd $0xff, %xmm11, %xmm11

  pshufd $0x00, %xmm15, %xmm12
  pshufd $0x55, %xmm15, %xmm13
  pshufd $0xaa, %xmm15, %xmm14
  pshufd $0xff, %xmm15, %xmm15

  /* Copy reordered A elements into buffer. */
  movaps %xmm0,   0*4*4(%rsp) /* A_{11} */
  movaps %xmm1,   1*4*4(%rsp) /* A_{12} */
  movaps %xmm2,   2*4*4(%rsp) /* A_{13} */
  movaps %xmm3,   3*4*4(%rsp) /* A_{14} */

  movaps %xmm4,   4*4*4(%rsp) /* A_{21} */
  movaps %xmm5,   5*4*4(%rsp) /* A_{22} */
  movaps %xmm6,   6*4*4(%rsp) /* A_{23} */
  movaps %xmm7,   7*4*4(%rsp) /* A_{24} */

  movaps %xmm8,   8*4*4(%rsp) /* A_{31} */
  movaps %xmm9,   9*4*4(%rsp) /* A_{32} */
  movaps %xmm10, 10*4*4(%rsp) /* A_{33} */
  movaps %xmm11, 11*4*4(%rsp) /* A_{34} */

  movaps %xmm12, 12*4*4(%rsp) /* A_{41} */
  movaps %xmm13, 13*4*4(%rsp) /* A_{42} */
  movaps %xmm14, 14*4*4(%rsp) /* A_{43} */
  movaps %xmm15, 15*4*4(%rsp) /* A_{44} */

  /* Zero out result matrix. */
  pxor %xmm8,  %xmm8
  pxor %xmm9,  %xmm9
  pxor %xmm10, %xmm10
  pxor %xmm11, %xmm11

  /* Load elements. */
  movaps  0*4*4(B),    %xmm0 /* B_{11} B_{12} B_{13} B_{14} */
  movaps  0*4*4(%rsp), %xmm2 /* A_{11} A_{11} A_{11} A_{11} */
  movaps  4*4*4(%rsp), %xmm3 /* A_{21} A_{21} A_{21} A_{21} */
  movaps  8*4*4(%rsp), %xmm4 /* A_{31} A_{31} A_{31} A_{31} */
  movaps 12*4*4(%rsp), %xmm5 /* A_{41} A_{41} A_{41} A_{41} */

  /* Multiply A*B and store in C. */
  mulps %xmm0, %xmm2
  addps %xmm2, %xmm8
  mulps %xmm0, %xmm3
  addps %xmm3, %xmm9
  mulps %xmm0, %xmm4
  addps %xmm4, %xmm10
  mulps %xmm0, %xmm5
  addps %xmm5, %xmm11

  /* Load elements. */
  movaps  1*4*4(B),    %xmm0 /* B_{21} B_{22} B_{23} B_{24} */
  movaps  1*4*4(%rsp), %xmm2 /* A_{12} A_{12} A_{12} A_{12} */
  movaps  5*4*4(%rsp), %xmm3 /* A_{22} A_{22} A_{22} A_{22} */
  movaps  9*4*4(%rsp), %xmm4 /* A_{32} A_{32} A_{32} A_{32} */
  movaps 13*4*4(%rsp), %xmm5 /* A_{42} A_{42} A_{42} A_{42} */

  /* Multiply A*B and store in C. */
  mulps %xmm0, %xmm2
  addps %xmm2, %xmm8
  mulps %xmm0, %xmm3
  addps %xmm3, %xmm9
  mulps %xmm0, %xmm4
  addps %xmm4, %xmm10
  mulps %xmm0, %xmm5
  addps %xmm5, %xmm11

  /* Load elements. */
  movaps  2*4*4(B),    %xmm0 /* B_{31} B_{32} B_{33} B_{34} */
  movaps  2*4*4(%rsp), %xmm2 /* A_{13} A_{13} A_{13} A_{13} */
  movaps  6*4*4(%rsp), %xmm3 /* A_{23} A_{23} A_{23} A_{23} */
  movaps 10*4*4(%rsp), %xmm4 /* A_{33} A_{33} A_{33} A_{33} */
  movaps 14*4*4(%rsp), %xmm5 /* A_{43} A_{43} A_{43} A_{43} */

  /* Multiply A*B and store in C. */
  mulps %xmm0, %xmm2
  addps %xmm2, %xmm8
  mulps %xmm0, %xmm3
  addps %xmm3, %xmm9
  mulps %xmm0, %xmm4
  addps %xmm4, %xmm10
  mulps %xmm0, %xmm5
  addps %xmm5, %xmm11

  /* Load elements. */
  movaps  3*4*4(B),    %xmm0 /* B_{41} B_{42} B_{43} B_{44} */
  movaps  3*4*4(%rsp), %xmm2 /* A_{14} A_{14} A_{14} A_{14} */
  movaps  7*4*4(%rsp), %xmm3 /* A_{24} A_{24} A_{24} A_{24} */
  movaps 11*4*4(%rsp), %xmm4 /* A_{34} A_{34} A_{34} A_{34} */
  movaps 15*4*4(%rsp), %xmm5 /* A_{44} A_{44} A_{44} A_{44} */

  /* Multiply A*B and store in C. */
  mulps %xmm0, %xmm2
  addps %xmm2, %xmm8
  mulps %xmm0, %xmm3
  addps %xmm3, %xmm9
  mulps %xmm0, %xmm4
  addps %xmm4, %xmm10
  mulps %xmm0, %xmm5
  addps %xmm5, %xmm11

  /* Add C from registers back into memory. */
  movaps 0*4*4(C), %xmm0
  movaps 1*4*4(C), %xmm1
  movaps 2*4*4(C), %xmm2
  movaps 3*4*4(C), %xmm3

  addps %xmm0, %xmm8
  addps %xmm1, %xmm9
  addps %xmm2, %xmm10
  addps %xmm3, %xmm11

  movaps %xmm8,  0*4*4(C)
  movaps %xmm9,  1*4*4(C)
  movaps %xmm10, 2*4*4(C)
  movaps %xmm11, 3*4*4(C)

  /* Restore stack. */
  movq OLD_STACK, %rsp
  ret

  .size sgemm_kernel_1, .-sgemm_kernel_1
