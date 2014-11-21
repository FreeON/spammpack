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
#define mask          %r13

#define stacksize 0x80
#define alpha     (%rsp)
#define old_MXCSR 0x40(%rsp)
#define new_MXCSR 0x44(%rsp)

#define SIZE 8 /* sizeof(void*) */

  .text
  .global stream_kernel_13
  .type stream_kernel_13, @function

  /* Function declaration in C.
   *
   *  struct multiply_stream_t
   *  {
   *    float *A_block;
   *    float *B_block;
   *    float *C_block;
   *    char  mask[8];
   *  };
   *
   *  void
   *  stream_kernel_13 (const unsigned int number_stream_elements,
   *      float alpha,
   *      struct multiply_stream_t *multiply_stream);
   */

  .align 256
stream_kernel_13:

  /* Save local variables on stack. */
  push old_stack
  push temp_reg_64
  push index
  push index_pointer
  push A
  push B
  push C
  push mask

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

  /* Divide number_stream_elements by 8 (to get a stride of 8 and to avoid
   * having to stride through the stream which incurs a higher latency. A
   * simple stride in the stream is faster.). */
  sar $3, number_stream_elements

  .align 16
loop:

  /* Update pointer to A. */
  mov (index_pointer, multiply_stream, 1), A
  mov 0x8(index_pointer, multiply_stream, 1), B
  mov 0x10(index_pointer, multiply_stream, 1), C
  lea 0x18(index_pointer, multiply_stream, 1), mask

  /* Increment loop index. */
  add $1, index
  //add $32, index_pointer

  /* Multiply block A(1,1)*B(1,1). */
  cmpb $0, 0(mask)
  je multiply_block_1

  movaps  0*4*4(A), %xmm0  /* A(1,1) */
  movaps  1*4*4(A), %xmm1  /* A(1,2) */
  movaps  2*4*4(A), %xmm2  /* A(1,3) */
  movaps  3*4*4(A), %xmm3  /* A(1,4) */
  movaps  0*4*4(B), %xmm4  /* B(1,1:4) */

  movaps  4*4*4(A), %xmm5  /* A(2,1) */
  movaps  5*4*4(A), %xmm6  /* A(2,2) */
  movaps  6*4*4(A), %xmm7  /* A(2,3) */
  movaps  7*4*4(A), %xmm8  /* A(2,4) */
  movaps  1*4*4(B), %xmm9  /* B(2,1:4) */

  movaps  8*4*4(A), %xmm10 /* A(3,1) */
  movaps  9*4*4(A), %xmm11 /* A(3,2) */
  movaps 10*4*4(A), %xmm12 /* A(3,3) */
  movaps 11*4*4(A), %xmm13 /* A(3,4) */
  movaps  2*4*4(B), %xmm14 /* B(3,1:4) */

  movaps 12*4*4(A), %xmm15 /* A(4,1) */
  movaps 13*4*4(A), %xmm0  /* A(4,2) */
  movaps 14*4*4(A), %xmm1  /* A(4,3) */
  movaps 15*4*4(A), %xmm2  /* A(4,4) */
  movaps  3*4*4(B), %xmm3  /* B(4,1:4) */

  .align 16
multiply_block_1:

  /* Multiply block A(1,2)*B(2,1). */
  cmpb $0, 1(mask)
  je multiply_block_2

  movaps 256+0*4*4(A), %xmm5  /* A(1,1) */
  movaps 256+1*4*4(A), %xmm6  /* A(1,2) */
  movaps 256+2*4*4(A), %xmm7  /* A(1,3) */
  movaps 256+3*4*4(A), %xmm8  /* A(1,4) */
  movaps 128+0*4*4(B), %xmm9  /* B(1,1:4) */

  movaps 256+4*4*4(A), %xmm10 /* A(2,1) */
  movaps 256+5*4*4(A), %xmm11 /* A(2,2) */
  movaps 256+6*4*4(A), %xmm12 /* A(2,3) */
  movaps 256+7*4*4(A), %xmm13 /* A(2,4) */
  movaps 128+1*4*4(B), %xmm14 /* B(2,1:4) */

  movaps 256+ 8*4*4(A), %xmm15 /* A(3,1) */
  movaps 256+ 9*4*4(A), %xmm0  /* A(3,2) */
  movaps 256+10*4*4(A), %xmm1  /* A(3,3) */
  movaps 256+11*4*4(A), %xmm2  /* A(3,4) */
  movaps 128+ 2*4*4(B), %xmm3  /* B(3,1:4) */

  movaps 256+12*4*4(A), %xmm4  /* A(4,1) */
  movaps 256+13*4*4(A), %xmm5  /* A(4,2) */
  movaps 256+14*4*4(A), %xmm6  /* A(4,3) */
  movaps 256+15*4*4(A), %xmm7  /* A(4,4) */
  movaps 128+ 3*4*4(B), %xmm8  /* B(4,1:4) */

  movaps 0*4*4(C), %xmm9 /* C(1,1:4) */

  .align 16
multiply_block_2:

  /* Multiply block A(1,1)*B(1,2). */
  cmpb $0, 2(mask)
  je multiply_block_3

  movaps  0*4*4(A), %xmm10 /* A(1,1) */
  movaps  1*4*4(A), %xmm11 /* A(1,2) */
  movaps  2*4*4(A), %xmm12 /* A(1,3) */
  movaps  3*4*4(A), %xmm13 /* A(1,4) */
  movaps  64+0*4*4(B), %xmm13 /* B(1,1:4) */

  movaps  4*4*4(A), %xmm15 /* A(2,1) */
  movaps  5*4*4(A), %xmm0  /* A(2,2) */
  movaps  6*4*4(A), %xmm1  /* A(2,3) */
  movaps  7*4*4(A), %xmm2  /* A(2,4) */
  movaps  64+1*4*4(B), %xmm3  /* B(2,1:4) */

  movaps  8*4*4(A), %xmm4  /* A(3,1) */
  movaps  9*4*4(A), %xmm5  /* A(3,2) */
  movaps 10*4*4(A), %xmm6  /* A(3,3) */
  movaps 11*4*4(A), %xmm7  /* A(3,4) */
  movaps  64+2*4*4(B), %xmm8  /* B(3,1:4) */

  movaps 12*4*4(A), %xmm9  /* A(4,1) */
  movaps 13*4*4(A), %xmm10 /* A(4,2) */
  movaps 14*4*4(A), %xmm11 /* A(4,3) */
  movaps 15*4*4(A), %xmm12 /* A(4,4) */
  movaps  64+3*4*4(B), %xmm13 /* B(4,1:4) */

  .align 16
multiply_block_3:

  /* Multiply block A(1,2)*B(2,2). */
  cmpb $0, 3(mask)
  je multiply_block_4

  movaps 256+0*4*4(A), %xmm14 /* A(1,1) */
  movaps 256+1*4*4(A), %xmm15 /* A(1,2) */
  movaps 256+2*4*4(A), %xmm0  /* A(1,3) */
  movaps 256+3*4*4(A), %xmm1  /* A(1,4) */
  movaps 192+0*4*4(B), %xmm2  /* B(1,1:4) */

  movaps 256+4*4*4(A), %xmm3  /* A(2,1) */
  movaps 256+5*4*4(A), %xmm4  /* A(2,2) */
  movaps 256+6*4*4(A), %xmm5  /* A(2,3) */
  movaps 256+7*4*4(A), %xmm6  /* A(2,4) */
  movaps 192+1*4*4(B), %xmm7  /* B(2,1:4) */

  movaps 256+ 8*4*4(A), %xmm8  /* A(3,1) */
  movaps 256+ 9*4*4(A), %xmm9  /* A(3,2) */
  movaps 256+10*4*4(A), %xmm10 /* A(3,3) */
  movaps 256+11*4*4(A), %xmm11 /* A(3,4) */
  movaps 192+ 2*4*4(B), %xmm12 /* B(3,1:4) */

  movaps 256+12*4*4(A), %xmm13 /* A(4,1) */
  movaps 256+13*4*4(A), %xmm14 /* A(4,2) */
  movaps 256+14*4*4(A), %xmm15 /* A(4,3) */
  movaps 256+15*4*4(A), %xmm0  /* A(4,4) */
  movaps 192+ 3*4*4(B), %xmm1  /* B(4,1:4) */

  movaps 1*4*4(C), %xmm2  /* C(2,1:4) */

  .align 16
multiply_block_4:

  /* Multiply block A(2,1)*B(1,1). */
  cmpb $0, 4(mask)
  je multiply_block_5

  movaps 512+ 0*4*4(A), %xmm3  /* A(1,1) */
  movaps 512+ 1*4*4(A), %xmm4  /* A(1,2) */
  movaps 512+ 2*4*4(A), %xmm5  /* A(1,3) */
  movaps 512+ 3*4*4(A), %xmm6  /* A(1,4) */
  movaps 0*4*4(B), %xmm7  /* B(1,1:4) */

  movaps 512+ 4*4*4(A), %xmm8  /* A(2,1) */
  movaps 512+ 5*4*4(A), %xmm9  /* A(2,2) */
  movaps 512+ 6*4*4(A), %xmm10 /* A(2,3) */
  movaps 512+ 7*4*4(A), %xmm11 /* A(2,4) */
  movaps 1*4*4(B), %xmm12 /* B(2,1:4) */

  movaps 512+ 8*4*4(A), %xmm13 /* A(3,1) */
  movaps 512+ 9*4*4(A), %xmm14 /* A(3,2) */
  movaps 512+10*4*4(A), %xmm15 /* A(3,3) */
  movaps 512+11*4*4(A), %xmm0  /* A(3,4) */
  movaps 2*4*4(B), %xmm1  /* B(3,1:4) */

  movaps 512+12*4*4(A), %xmm2  /* A(4,1) */
  movaps 512+13*4*4(A), %xmm3  /* A(4,2) */
  movaps 512+14*4*4(A), %xmm4  /* A(4,3) */
  movaps 512+15*4*4(A), %xmm5  /* A(4,4) */
  movaps 3*4*4(B), %xmm6  /* B(4,1:4) */

  .align 16
multiply_block_5:

  /* Multiply block A(2,2)*B(2,1). */
  cmpb $0, 5(mask)
  je multiply_block_6

  movaps 768+ 0*4*4(A), %xmm7  /* A(1,1) */
  movaps 768+ 1*4*4(A), %xmm8  /* A(1,2) */
  movaps 768+ 2*4*4(A), %xmm9  /* A(1,3) */
  movaps 768+ 3*4*4(A), %xmm10 /* A(1,4) */
  movaps 128+0*4*4(B), %xmm11 /* B(1,1:4) */

  movaps 768+ 4*4*4(A), %xmm12 /* A(2,1) */
  movaps 768+ 5*4*4(A), %xmm13 /* A(2,2) */
  movaps 768+ 6*4*4(A), %xmm14 /* A(2,3) */
  movaps 768+ 7*4*4(A), %xmm15 /* A(2,4) */
  movaps 128+1*4*4(B), %xmm0  /* B(2,1:4) */

  movaps 768+ 8*4*4(A), %xmm1  /* A(3,1) */
  movaps 768+ 9*4*4(A), %xmm2  /* A(3,2) */
  movaps 768+10*4*4(A), %xmm3  /* A(3,3) */
  movaps 768+11*4*4(A), %xmm4  /* A(3,4) */
  movaps 128+2*4*4(B), %xmm5  /* B(3,1:4) */

  movaps 768+12*4*4(A), %xmm6  /* A(4,1) */
  movaps 768+13*4*4(A), %xmm7  /* A(4,2) */
  movaps 768+14*4*4(A), %xmm8  /* A(4,3) */
  movaps 768+15*4*4(A), %xmm9  /* A(4,4) */
  movaps 128+3*4*4(B), %xmm10 /* B(4,1:4) */

  movaps 2*4*4(C), %xmm11 /* C(3,1:4) */

  .align 16
multiply_block_6:

  /* Multiply block A(2,1)*B(1,2). */
  cmpb $0, 6(mask)
  je multiply_block_7

  movaps 512+ 0*4*4(A), %xmm12 /* A(1,1) */
  movaps 512+ 1*4*4(A), %xmm13 /* A(1,2) */
  movaps 512+ 2*4*4(A), %xmm14 /* A(1,3) */
  movaps 512+ 3*4*4(A), %xmm15 /* A(1,4) */
  movaps  64+0*4*4(B), %xmm0  /* B(1,1:4) */

  movaps 512+ 4*4*4(A), %xmm1  /* A(2,1) */
  movaps 512+ 5*4*4(A), %xmm2  /* A(2,2) */
  movaps 512+ 6*4*4(A), %xmm3  /* A(2,3) */
  movaps 512+ 7*4*4(A), %xmm4  /* A(2,4) */
  movaps  64+1*4*4(B), %xmm5  /* B(2,1:4) */

  movaps 512+ 8*4*4(A), %xmm6  /* A(3,1) */
  movaps 512+ 9*4*4(A), %xmm7  /* A(3,2) */
  movaps 512+10*4*4(A), %xmm8  /* A(3,3) */
  movaps 512+11*4*4(A), %xmm9  /* A(3,4) */
  movaps  64+2*4*4(B), %xmm10 /* B(3,1:4) */

  movaps 512+12*4*4(A), %xmm11 /* A(4,1) */
  movaps 512+13*4*4(A), %xmm12 /* A(4,2) */
  movaps 512+14*4*4(A), %xmm13 /* A(4,3) */
  movaps 512+15*4*4(A), %xmm14 /* A(4,4) */
  movaps  64+3*4*4(B), %xmm15 /* B(4,1:4) */

  .align 16
multiply_block_7:

  /* Multiply block A(2,2)*B(2,2). */
  cmpb $0, 7(mask)
  je done_multiply

  movaps 768+ 0*4*4(A), %xmm0  /* A(1,1) */
  movaps 768+ 1*4*4(A), %xmm1  /* A(1,2) */
  movaps 768+ 2*4*4(A), %xmm2  /* A(1,3) */
  movaps 768+ 3*4*4(A), %xmm3  /* A(1,4) */
  movaps 192+0*4*4(B), %xmm4  /* B(1,1:4) */

  movaps 768+ 4*4*4(A), %xmm5  /* A(2,1) */
  movaps 768+ 5*4*4(A), %xmm6  /* A(2,2) */
  movaps 768+ 6*4*4(A), %xmm7  /* A(2,3) */
  movaps 768+ 7*4*4(A), %xmm8  /* A(2,4) */
  movaps 192+1*4*4(B), %xmm9  /* B(2,1:4) */

  movaps 768+ 8*4*4(A), %xmm10 /* A(3,1) */
  movaps 768+ 9*4*4(A), %xmm11 /* A(3,2) */
  movaps 768+10*4*4(A), %xmm12 /* A(3,3) */
  movaps 768+11*4*4(A), %xmm13 /* A(3,4) */
  movaps 192+2*4*4(B), %xmm14 /* B(3,1:4) */

  movaps 768+12*4*4(A), %xmm15 /* A(4,1) */
  movaps 768+13*4*4(A), %xmm0  /* A(4,2) */
  movaps 768+14*4*4(A), %xmm1  /* A(4,3) */
  movaps 768+15*4*4(A), %xmm2  /* A(4,4) */
  movaps 192+3*4*4(B), %xmm3  /* B(4,1:4) */

  movaps 3*4*4(C), %xmm4  /* C(4,1:4) */

  .align 16
done_multiply:

  /* Loop comparison. */
  cmp number_stream_elements, index

  /* Finish loop. */
  jl loop

  .align 16
done:

  /* Restore rounding mode. */
  ldmxcsr old_MXCSR

  /* Restore old stack. */
  mov old_stack, %rsp

  /* Restore local variables from stack. */
  pop mask
  pop C
  pop B
  pop A
  pop index_pointer
  pop index
  pop temp_reg_64
  pop old_stack

  ret

  .size stream_kernel_13, .-stream_kernel_13
