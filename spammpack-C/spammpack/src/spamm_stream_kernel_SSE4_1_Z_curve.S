# Function ABI.
#define number_stream_elements %rdi
#define alpha                  %xmm0
#define tolerance              %xmm1
#define multiply_stream        %rsi

# Define loop variables.
#define index        %rax
#define base_pointer %rdx

# Define a variable for the old stack.
#define old_stack %rcx

# Define jump table index.
#define jump_index %r11
#define jump_index_temp %r12

# Define pointers to matrix data nodes in stream.
#define A %r8
#define B %r9
#define C %r10

  # Function prolog.
  .text
  .align 256
  .global spamm_stream_kernel_SSE4_1_Z_curve
  .type spamm_stream_kernel_SSE4_1_Z_curve, @function

spamm_stream_kernel_SSE4_1_Z_curve:

  # Push used registers on stack.
  push index
  push base_pointer
  push A
  push B
  push C
  push jump_index
  push jump_index_temp

  # Push stack pointer so we can make room for local storage.
  push old_stack
  mov %rsp, old_stack
  sub $0x16, %rsp
  and $-0x10, %rsp

  # Copy alpha into all 4 elements of SSE register.
  shufps $0x0, alpha, alpha

  # Test whether number_stream_elements is zero.
  test number_stream_elements, number_stream_elements
  jbe stream_done

  # Set loop index to zero.
  xor base_pointer, base_pointer
  xor index, index

  .align 16
stream_loop:

  # Set the base pointer using sizeof(multiply_stream_t) = 24 (0x18).
  imul $0x18, base_pointer, base_pointer

  # Load pointers to stream matrix blocks.
  mov (multiply_stream, base_pointer, 1), A
  mov 0x8(multiply_stream, base_pointer, 1), B
  mov 0x10(multiply_stream, base_pointer, 1), C

  # Load shuffle mask for norm comparisons.
  movl $0x0c080400, 0x0(%rsp)
  movl $0x80808080, 0x4(%rsp)
  movl $0x80808080, 0x8(%rsp)
  movl $0x80808080, 0xc(%rsp)
  movaps (%rsp), %xmm2

  # First level of hierarchy. Do some norm products and decide where to go on
  # the next tier.
  movaps 0x10(A), %xmm3
  mulps 0x10(B), %xmm3
  cmpps $0x02, tolerance, %xmm3 # $xmm3 <= tolerance?
  pshufb %xmm2, %xmm3
  pmovmskb %xmm3, jump_index

  # Jump table for first tier.
  lea (,jump_index, 4), jump_index
  lea jump_table(%rip), jump_index_temp
  mov (jump_index_temp, jump_index), jump_index
  lea jump_table(%rip), jump_index_temp
  lea (jump_index_temp, jump_index), jump_index
  jmp *jump_index

  .section .rodata
  .align 4
jump_table:
  .long tier_00-jump_table
  .long tier_01-jump_table
  .long tier_02-jump_table
  .long tier_03-jump_table
  .long tier_04-jump_table
  .long tier_05-jump_table
  .long tier_06-jump_table
  .long tier_07-jump_table
  .long tier_08-jump_table

  .text
tier_00:
  jmp loop_end
tier_01:
  jmp loop_end
tier_02:
  jmp loop_end
tier_03:
  jmp loop_end
tier_04:
  jmp loop_end
tier_05:
  jmp loop_end
tier_06:
  jmp loop_end
tier_07:
  jmp loop_end
tier_08:
  jmp loop_end

loop_end:
  # Loop end.
  add $0x01, index
  mov index, base_pointer
  cmp number_stream_elements, index
  jb stream_loop

  .align 16
stream_done:

  # Pop registers from stack.
  mov old_stack, %rsp
  pop old_stack
  pop jump_index_temp
  pop jump_index
  pop C
  pop B
  pop A
  pop base_pointer
  pop index

  # Return from function.
  ret

  # Function epilog.
  .size spamm_stream_kernel_SSE4_1_Z_curve, .-spamm_stream_kernel_SSE4_1_Z_curve
