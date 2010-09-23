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
#define A             %r10
#define B             %r11
#define C             %r12
#define A_ahead       %r13
#define B_ahead       %r14
#define C_ahead       %r15

#define stacksize 0x80
#define alpha     (%rsp)
#define old_MXCSR 0x40(%rsp)
#define new_MXCSR 0x44(%rsp)

#define SIZE 8 /* sizeof(void*) */

  .text
  .global stream_kernel_5
  .type stream_kernel_5, @function

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
   *  stream_kernel_5 (const unsigned int number_stream_elements,
   *      float alpha,
   *      struct multiply_stream_t *multiply_stream);
   */

  .align 256
stream_kernel_5:

  /* Save local variables on stack. */
  push old_stack
  push temp_reg_64
  push index
  push index_pointer
  push A
  push B
  push C
  push A_ahead
  push B_ahead
  push C_ahead

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

  .align 16
loop:

  mov (index_pointer, multiply_stream, 1), A
  mov 0x8(index_pointer, multiply_stream, 1), B
  mov 0x10(index_pointer, multiply_stream, 1), C

  /* Increment loop index. */
  add $1, index
  add $0x18, index_pointer

  /* Loop comparison. */
  cmp number_stream_elements, index

  /* Multiply block. */

  /* Do nothing. */

  /* Finish loop. */
  jl loop

  .align 16
done:

  /* Restore rounding mode. */
  ldmxcsr old_MXCSR

  /* Restore old stack. */
  mov old_stack, %rsp

  /* Restore local variables from stack. */
  pop C_ahead
  pop B_ahead
  pop A_ahead
  pop C
  pop B
  pop A
  pop index_pointer
  pop index
  pop temp_reg_64
  pop old_stack

  ret

  .size stream_kernel_5, .-stream_kernel_5
