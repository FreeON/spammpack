/* Define function arguments. */
#define number_stream_elements  %rdi
#define alpha_arg               %xmm0
#define multiply_stream         %rsi

/* Define local variables. */
#define old_stack     %r8
#define index_outer   %r9
#define index_inner   %r10
#define A             %r11
#define B             %r12
#define C             %r13
#define temp_reg_64   %rcx
#define temp_reg_32   %ecx

#define stacksize 0x80
#define alpha     (%rsp)
#define old_MXCSR 0x40(%rsp)
#define new_MXCSR 0x44(%rsp)

#define SIZE 8 /* sizeof(void*) */

  .text
  .global stream_kernel_9
  .type stream_kernel_9, @function

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
   *  stream_kernel_9 (const unsigned int number_stream_elements,
   *      float alpha,
   *      struct multiply_stream_t *multiply_stream);
   */

  .align 256
stream_kernel_9:

  /* Save local variables on stack. */
  push old_stack
  push temp_reg_64
  push index_inner
  push index_outer
  push A
  push B
  push C
  push %rax
  push %rdx

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

  /* index_outer = 0 */
  xor index_outer, index_outer

  .align 16
loop_outer:

  /* Load pointer to beginning of stream. */
  mov (multiply_stream), A

  /* Increment loop index. */
  add $16, index_outer

  /* index_inner = 0 */
  xor index_inner, index_inner

  .align 16
loop_inner:

  /* Increment loop index. */
  add $1, index_inner

  /* Multiply block. */
  movaps  0*4*4(A), %xmm0  /* A(1,1) */
  movaps  1*4*4(A), %xmm1  /* A(1,2) */
  movaps  2*4*4(A), %xmm2  /* A(1,3) */
  movaps  3*4*4(A), %xmm3  /* A(1,4) */
  movaps  4*4*4(A), %xmm4  /* A(2,1) */
  movaps  5*4*4(A), %xmm5  /* A(2,2) */
  movaps  6*4*4(A), %xmm6  /* A(2,3) */
  movaps  7*4*4(A), %xmm7  /* A(2,4) */
  movaps  8*4*4(A), %xmm8  /* A(3,1) */
  movaps  9*4*4(A), %xmm9  /* A(3,2) */
  movaps 10*4*4(A), %xmm10 /* A(3,3) */
  movaps 11*4*4(A), %xmm11 /* A(3,4) */
  movaps 12*4*4(A), %xmm12 /* A(4,1) */
  movaps 13*4*4(A), %xmm13 /* A(4,2) */
  movaps 14*4*4(A), %xmm14 /* A(4,3) */
  movaps 15*4*4(A), %xmm15 /* A(4,4) */

  /* Increment A pointer to next block. */
  add $256, A

  /* Loop comparison. */
  cmp $16, index_inner

  /* Finish loop. */
  jl loop_inner

  /* Loop comparison. */
  cmp number_stream_elements, index_outer

  /* Finish loop. */
  jl loop_outer

  .align 16
done:

  /* Restore rounding mode. */
  ldmxcsr old_MXCSR

  /* Restore old stack. */
  mov old_stack, %rsp

  /* Restore local variables from stack. */
  pop %rdx
  pop %rax
  pop C
  pop B
  pop A
  pop index_outer
  pop index_inner
  pop temp_reg_64
  pop old_stack

  ret

  .size stream_kernel_9, .-stream_kernel_9
