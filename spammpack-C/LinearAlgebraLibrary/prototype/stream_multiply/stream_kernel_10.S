/* Define function arguments. */
#define number_stream_elements  %rdi
#define alpha_arg               %xmm0
#define multiply_stream         %rsi

/* Define local variables. */
#define old_stack     %r8
#define index         %rax
#define index_pointer %rdx
#define temp_reg_64   %rcx
#define temp_reg_32   %ecx
#define A_0           %r10
#define A_1           %r11
#define A_2           %r12
#define A_3           %r13

#define stacksize 0x80
#define alpha     (%rsp)
#define old_MXCSR 0x40(%rsp)
#define new_MXCSR 0x44(%rsp)

#define SIZE 8 /* sizeof(void*) */

  .text
  .global stream_kernel_10
  .type stream_kernel_10, @function

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
   *  stream_kernel_10 (const unsigned int number_stream_elements,
   *      float alpha,
   *      struct multiply_stream_t *multiply_stream);
   */

  .align 256
stream_kernel_10:

  /* Save local variables on stack. */
  push old_stack
  push temp_reg_64
  push index
  push index_pointer
  push A_0
  push A_1
  push A_2
  push A_3

  mov %rsp, old_stack /* Save old stack. */
  sub $stacksize, %rsp /* Make some room for local variables. */
  and $-0x40, %rsp /* Align stack on 64-byte boundary (cache line). */

  /* Do we have to loop at all? */
  test number_stream_elements, number_stream_elements
  je done

  /* Dilate alpha to full SSE vector. */
  pshufd $0x00, alpha_arg, alpha_arg
  movaps alpha_arg, alpha

  /* Set rounding mode. */
  stmxcsr old_MXCSR
  mov old_MXCSR, temp_reg_32
  or $0x8040, temp_reg_32
  mov temp_reg_32, new_MXCSR
  ldmxcsr new_MXCSR

  /* index = 0 */
  xor index, index
  xor index_pointer, index_pointer

  mov (index_pointer, multiply_stream, 1), A_0
  mov 24(index_pointer, multiply_stream, 1), A_1
  mov 48(index_pointer, multiply_stream, 1), A_2
  mov 72(index_pointer, multiply_stream, 1), A_3

  .align 16
loop:

  /* Increment loop index. */
  add $4, index
  add $96, index_pointer

  /* Loop comparison. */
  cmp number_stream_elements, index

  /* Multiply block. */
  movaps  0*4*4(A_0), %xmm0  /* A(1,1) */
  movaps  1*4*4(A_0), %xmm1  /* A(1,2) */
  movaps  2*4*4(A_0), %xmm2  /* A(1,3) */
  movaps  3*4*4(A_0), %xmm3  /* A(1,4) */
  movaps  4*4*4(A_0), %xmm4  /* A(2,1) */
  movaps  5*4*4(A_0), %xmm5  /* A(2,2) */
  movaps  6*4*4(A_0), %xmm6  /* A(2,3) */
  movaps  7*4*4(A_0), %xmm7  /* A(2,4) */
  movaps  8*4*4(A_0), %xmm8  /* A(3,1) */
  movaps  9*4*4(A_0), %xmm9  /* A(3,2) */
  movaps 10*4*4(A_0), %xmm10 /* A(3,3) */
  movaps 11*4*4(A_0), %xmm11 /* A(3,4) */
  movaps 12*4*4(A_0), %xmm12 /* A(4,1) */
  movaps 13*4*4(A_0), %xmm13 /* A(4,2) */
  movaps 14*4*4(A_0), %xmm14 /* A(4,3) */
  movaps 15*4*4(A_0), %xmm15 /* A(4,4) */

  mov (index_pointer, multiply_stream, 1), A_0

  /* Multiply block. */
  movaps  0*4*4(A_1), %xmm0  /* A(1,1) */
  movaps  1*4*4(A_1), %xmm1  /* A(1,2) */
  movaps  2*4*4(A_1), %xmm2  /* A(1,3) */
  movaps  3*4*4(A_1), %xmm3  /* A(1,4) */
  movaps  4*4*4(A_1), %xmm4  /* A(2,1) */
  movaps  5*4*4(A_1), %xmm5  /* A(2,2) */
  movaps  6*4*4(A_1), %xmm6  /* A(2,3) */
  movaps  7*4*4(A_1), %xmm7  /* A(2,4) */
  movaps  8*4*4(A_1), %xmm8  /* A(3,1) */
  movaps  9*4*4(A_1), %xmm9  /* A(3,2) */
  movaps 10*4*4(A_1), %xmm10 /* A(3,3) */
  movaps 11*4*4(A_1), %xmm11 /* A(3,4) */
  movaps 12*4*4(A_1), %xmm12 /* A(4,1) */
  movaps 13*4*4(A_1), %xmm13 /* A(4,2) */
  movaps 14*4*4(A_1), %xmm14 /* A(4,3) */
  movaps 15*4*4(A_1), %xmm15 /* A(4,4) */

  mov 24(index_pointer, multiply_stream, 1), A_1

  /* Multiply block. */
  movaps  0*4*4(A_2), %xmm0  /* A(1,1) */
  movaps  1*4*4(A_2), %xmm1  /* A(1,2) */
  movaps  2*4*4(A_2), %xmm2  /* A(1,3) */
  movaps  3*4*4(A_2), %xmm3  /* A(1,4) */
  movaps  4*4*4(A_2), %xmm4  /* A(2,1) */
  movaps  5*4*4(A_2), %xmm5  /* A(2,2) */
  movaps  6*4*4(A_2), %xmm6  /* A(2,3) */
  movaps  7*4*4(A_2), %xmm7  /* A(2,4) */
  movaps  8*4*4(A_2), %xmm8  /* A(3,1) */
  movaps  9*4*4(A_2), %xmm9  /* A(3,2) */
  movaps 10*4*4(A_2), %xmm10 /* A(3,3) */
  movaps 11*4*4(A_2), %xmm11 /* A(3,4) */
  movaps 12*4*4(A_2), %xmm12 /* A(4,1) */
  movaps 13*4*4(A_2), %xmm13 /* A(4,2) */
  movaps 14*4*4(A_2), %xmm14 /* A(4,3) */
  movaps 15*4*4(A_2), %xmm15 /* A(4,4) */

  mov 48(index_pointer, multiply_stream, 1), A_2

  /* Multiply block. */
  movaps  0*4*4(A_3), %xmm0  /* A(1,1) */
  movaps  1*4*4(A_3), %xmm1  /* A(1,2) */
  movaps  2*4*4(A_3), %xmm2  /* A(1,3) */
  movaps  3*4*4(A_3), %xmm3  /* A(1,4) */
  movaps  4*4*4(A_3), %xmm4  /* A(2,1) */
  movaps  5*4*4(A_3), %xmm5  /* A(2,2) */
  movaps  6*4*4(A_3), %xmm6  /* A(2,3) */
  movaps  7*4*4(A_3), %xmm7  /* A(2,4) */
  movaps  8*4*4(A_3), %xmm8  /* A(3,1) */
  movaps  9*4*4(A_3), %xmm9  /* A(3,2) */
  movaps 10*4*4(A_3), %xmm10 /* A(3,3) */
  movaps 11*4*4(A_3), %xmm11 /* A(3,4) */
  movaps 12*4*4(A_3), %xmm12 /* A(4,1) */
  movaps 13*4*4(A_3), %xmm13 /* A(4,2) */
  movaps 14*4*4(A_3), %xmm14 /* A(4,3) */
  movaps 15*4*4(A_3), %xmm15 /* A(4,4) */

  mov 72(index_pointer, multiply_stream, 1), A_3

  /* Finish loop. */
  jl loop

  .align 16
done:

  /* Restore rounding mode. */
  ldmxcsr old_MXCSR

  /* Restore old stack. */
  mov old_stack, %rsp

  /* Restore local variables from stack. */
  pop A_3
  pop A_2
  pop A_1
  pop A_0
  pop index_pointer
  pop index
  pop temp_reg_64
  pop old_stack

  ret

  .size stream_kernel_10, .-stream_kernel_10
