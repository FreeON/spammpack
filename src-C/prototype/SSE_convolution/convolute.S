  .text
  .global convolute
  .type convolute, @function

#define A_index_array %rdi
#define A_norm %rsi
#define B_index_array %rdx
#define B_norm %rcx

convolute:
  # Make some room on the stack.
  push %rax
  movq %rsp, %rax
  sub $0x10, %rsp
  and $-0x10, %rsp

  # Load one index value of A.
  movd 0(A_index_array), %xmm0
  shufps $0x00, %xmm0, %xmm0

  # Load first norm of A.
  movss 0(A_norm), %xmm1
  shufps $0x00, %xmm1, %xmm1

  # Load first 4 indices of B.
  movdqa 0(B_index_array), %xmm2

  # Restore stack.
  movq %rax, %rsp
  pop %rax

  # Return.
  ret

  .size convolute, .-convolute
