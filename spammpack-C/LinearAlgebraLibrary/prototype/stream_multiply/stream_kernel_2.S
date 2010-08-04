#define number_stream_elements  %rdi
#define alpha_arg               %xmm0
#define multiply_stream         %rsi

#define stream_element_size 3*8
#define buffer_size         4*4*4*4
#define stacksize           64

#define stream_index   %rax
#define stream_element %r8
#define old_stack      %rbx
#define alpha          buffer_size(%rsp)
#define A_dilated      %rsp
#define A              %r10
#define B              %r11
#define C              %r12

  .text
  .align 256
  .global stream_kernel_2
  .type stream_kernel_2, @function

#  Function declaration in C.
#
#  struct multiply_stream_t
#  {
#    float *A_block;
#    float *B_block;
#    float *C_block;
#  };
#
#  void
#  stream_kernel_2 (const unsigned int number_stream_elements,
#      float alpha,
#      struct multiply_stream_t *multiply_stream);

stream_kernel_2:

  /* Store clobbered registers on stack. */
  subq $stacksize, %rsp
  movq %rax,  0(%rsp)
  movq %r8,   8(%rsp)
  movq %rbx, 16(%rsp)
  movq %r10, 24(%rsp)
  movq %r11, 32(%rsp)
  movq %r12, 40(%rsp)
  movq %rsp, old_stack

  /* Create buffer on stack. */
  subq $(buffer_size+16), %rsp
  andq $-4096, %rsp /* Aling buffer on page boundary. */

  /* First we copy the first double word to the other double words in xmm0. */
  pshufd $0x00, alpha_arg, alpha_arg
  movaps alpha_arg, alpha

  /* Loop. */
  movq $0, stream_index
  leaq (multiply_stream), stream_element

loop:

  /* Calculate effective addresses of matrix blocks. */
  movq 0*8(stream_element), A
  leaq (A), A
  movq 1*8(stream_element), B
  leaq (B), B
  movq 2*8(stream_element), C
  leaq (C), C

  /* Dilate A into buffer. */
  movaps (A),     %xmm3
  movaps 4*4(A),  %xmm7
  movaps 8*4(A),  %xmm11
  movaps 12*4(A), %xmm15

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

  movaps %xmm0, 0*4*4(A_dilated)
  movaps %xmm1, 1*4*4(A_dilated)
  movaps %xmm2, 2*4*4(A_dilated)
  movaps %xmm3, 3*4*4(A_dilated)

  movaps %xmm4, 4*4*4(A_dilated)
  movaps %xmm5, 5*4*4(A_dilated)
  movaps %xmm6, 6*4*4(A_dilated)
  movaps %xmm7, 7*4*4(A_dilated)

  movaps %xmm8,   8*4*4(A_dilated)
  movaps %xmm9,   9*4*4(A_dilated)
  movaps %xmm10, 10*4*4(A_dilated)
  movaps %xmm11, 11*4*4(A_dilated)

  movaps %xmm12, 12*4*4(A_dilated)
  movaps %xmm13, 13*4*4(A_dilated)
  movaps %xmm14, 14*4*4(A_dilated)
  movaps %xmm15, 15*4*4(A_dilated)

  /* Multiply. */
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
  movaps alpha, %xmm0
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

  /* Increment pointer in multiply stream. */
  addq $3*8, stream_element

  /* Finish loop. */
  add $1, %rax
  cmp number_stream_elements, %rax
  jl loop

  /* We are done. */
  movq old_stack, %rsp
  movq 40(%rsp), %r12
  movq 32(%rsp), %r11
  movq 24(%rsp), %r10
  movq 16(%rsp), %rbx
  movq  8(%rsp), %r8
  movq  0(%rsp), %rax
  addq $stacksize, %rsp
  ret

  .size stream_kernel_2, .-stream_kernel_2
