#define number_stream_elements  %rdi
#define alpha_arg               %xmm0
#define multiply_stream         %rsi

#define stacksize  64

#define index     %rax
#define end_index %rcx
#define old_stack %rdx
#define A         %r10
#define B         %r11
#define C         %r12

#define A_dilated %r13

  .text
  .align 256
  .global stream_kernel_3
  .type stream_kernel_3, @function

  /* Function declaration in C.
   *
   *  struct multiply_stream_t
   *  {
   *    float *A_block;
   *    float *B_block;
   *    float *C_block;
   *  };
   *
   *  void
   *  stream_kernel_3 (const unsigned int number_stream_elements,
   *      float alpha,
   *      struct multiply_stream_t *multiply_stream);
   */

stream_kernel_3:

  subq $(stacksize), %rsp
  movq index,     0*8(%rsp)
  movq end_index, 1*8(%rsp)
  movq old_stack, 2*8(%rsp)
  movq A,         3*8(%rsp)
  movq B,         4*8(%rsp)
  movq C,         5*8(%rsp)
  movq A_dilated, 6*8(%rsp)

  movq %rsp, old_stack
  movq number_stream_elements, A_dilated
  shl $8, A_dilated
  lea 16( , A_dilated, 1), A_dilated
  subq A_dilated, %rsp
  andq $-4096, %rsp

  pshufd $0x00, alpha_arg, alpha_arg /* Copy the first double word to the other double words in xmm0. */
  movaps alpha_arg, (%rsp)

  test number_stream_elements, number_stream_elements
  je done

  xor index, index
  lea (number_stream_elements, number_stream_elements, 2), end_index
  shl $3, end_index

  lea 16(%rsp), A_dilated

  .align 16
loop_dilate:

  movq  0(multiply_stream, index, 1), A

  movaps (A),     %xmm3
  movaps 4*4(A),  %xmm7
  movaps 8*4(A),  %xmm11
  movaps 12*4(A), %xmm15

  pshufd $0x00, %xmm3,  %xmm0
  pshufd $0x55, %xmm3,  %xmm1
  pshufd $0xaa, %xmm3,  %xmm2
  pshufd $0xff, %xmm3,  %xmm3

  pshufd $0x00, %xmm7,  %xmm4
  pshufd $0x55, %xmm7,  %xmm5
  pshufd $0xaa, %xmm7,  %xmm6
  pshufd $0xff, %xmm7,  %xmm7

  pshufd $0x00, %xmm11, %xmm8
  pshufd $0x55, %xmm11, %xmm9
  pshufd $0xaa, %xmm11, %xmm10
  pshufd $0xff, %xmm11, %xmm11

  pshufd $0x00, %xmm15, %xmm12
  pshufd $0x55, %xmm15, %xmm13
  pshufd $0xaa, %xmm15, %xmm14
  pshufd $0xff, %xmm15, %xmm15

  movaps %xmm0,   0*4*4(A_dilated)
  movaps %xmm1,   1*4*4(A_dilated)
  movaps %xmm2,   2*4*4(A_dilated)
  movaps %xmm3,   3*4*4(A_dilated)

  movaps %xmm4,   4*4*4(A_dilated)
  movaps %xmm5,   5*4*4(A_dilated)
  movaps %xmm6,   6*4*4(A_dilated)
  movaps %xmm7,   7*4*4(A_dilated)

  movaps %xmm8,   8*4*4(A_dilated)
  movaps %xmm9,   9*4*4(A_dilated)
  movaps %xmm10, 10*4*4(A_dilated)
  movaps %xmm11, 11*4*4(A_dilated)

  movaps %xmm12, 12*4*4(A_dilated)
  movaps %xmm13, 13*4*4(A_dilated)
  movaps %xmm14, 14*4*4(A_dilated)
  movaps %xmm15, 15*4*4(A_dilated)

  addq $(256), A_dilated
  addq $(3*8), index
  cmp end_index, index
  jne loop_dilate

  xor index, index
  lea (number_stream_elements, number_stream_elements, 2), end_index
  shl $3, end_index

  lea 16(%rsp), A_dilated

  .align 16
loop_multiply:

  movq  8(multiply_stream, index, 1), B
  movq 16(multiply_stream, index, 1), C

  addq $(256), A_dilated
  addq $(3*8), index
  cmp end_index, index
  jne loop_multiply

  .align 16
done:

  movq old_stack, %rsp

  movq 6*8(%rsp), A_dilated
  movq 5*8(%rsp), C
  movq 4*8(%rsp), B
  movq 3*8(%rsp), A
  movq 2*8(%rsp), old_stack
  movq 1*8(%rsp), end_index
  movq 0*8(%rsp), index
  addq $(stacksize), %rsp

  ret

  .size stream_kernel_3, .-stream_kernel_3
