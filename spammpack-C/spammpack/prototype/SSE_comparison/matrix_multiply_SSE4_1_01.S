# C API:
#
# void
# matrix_multiply_SSE4_1_01 (const unsigned int N, float *A, float *B, float *C);

#define N %rdi
#define A %rsi
#define B %rdx
#define C %rcx

#define i %rax

  .text
  .align 256
  .global matrix_multiply_SSE4_1_01
  .type matrix_multiply_SSE4_1_01, @function

matrix_multiply_SSE4_1_01:

  push i
  xor i, i

  test N, N
  jbe end_loop

start_loop:

  movaps 0x00(C), %xmm0
  movaps 0x10(C), %xmm1
  movaps 0x20(C), %xmm2
  movaps 0x30(C), %xmm3

  movaps 0x00(B), %xmm4
  movaps 0x10(B), %xmm5
  movaps 0x20(B), %xmm6
  movaps 0x30(B), %xmm7

  movaps 0x00(A), %xmm8
  movaps 0x10(A), %xmm9

  # Calculate C(1,:).
  movaps %xmm4, %xmm10
  dpps $0xf1, %xmm8, %xmm10
  movaps %xmm5, %xmm11
  dpps $0xf2, %xmm8, %xmm11
  movaps %xmm6, %xmm12
  dpps $0xf4, %xmm8, %xmm12
  movaps %xmm7, %xmm13
  dpps $0xf8, %xmm8, %xmm13
  blendps $0x01, %xmm10, %xmm11
  blendps $0x03, %xmm11, %xmm12
  blendps $0x07, %xmm12, %xmm13
  addps %xmm13, %xmm0

  movaps 0x20(A), %xmm8

  # Calculate C(2,:).
  movaps %xmm4, %xmm10
  dpps $0xf1, %xmm9, %xmm10
  movaps %xmm5, %xmm11
  dpps $0xf2, %xmm9, %xmm11
  movaps %xmm6, %xmm12
  dpps $0xf4, %xmm9, %xmm12
  movaps %xmm7, %xmm13
  dpps $0xf8, %xmm9, %xmm13
  blendps $0x01, %xmm10, %xmm11
  blendps $0x03, %xmm11, %xmm12
  blendps $0x07, %xmm12, %xmm13
  addps %xmm13, %xmm1

  movaps 0x30(A), %xmm9

  # Calculate C(3,:).
  movaps %xmm4, %xmm10
  dpps $0xf1, %xmm8, %xmm10
  movaps %xmm5, %xmm11
  dpps $0xf2, %xmm8, %xmm11
  movaps %xmm6, %xmm12
  dpps $0xf4, %xmm8, %xmm12
  movaps %xmm7, %xmm13
  dpps $0xf8, %xmm8, %xmm13
  blendps $0x01, %xmm10, %xmm11
  blendps $0x03, %xmm11, %xmm12
  blendps $0x07, %xmm12, %xmm13
  addps %xmm13, %xmm2

  # Calculate C(4,:).
  movaps %xmm4, %xmm10
  dpps $0xf1, %xmm9, %xmm10
  movaps %xmm5, %xmm11
  dpps $0xf2, %xmm9, %xmm11
  movaps %xmm6, %xmm12
  dpps $0xf4, %xmm9, %xmm12
  movaps %xmm7, %xmm13
  dpps $0xf8, %xmm9, %xmm13
  blendps $0x01, %xmm10, %xmm11
  blendps $0x03, %xmm11, %xmm12
  blendps $0x07, %xmm12, %xmm13
  addps %xmm13, %xmm3

  # Write C back.
  movaps %xmm0, 0x00(C)
  movaps %xmm1, 0x10(C)
  movaps %xmm2, 0x20(C)
  movaps %xmm3, 0x30(C)

  inc i
  cmp N, i
  jb start_loop

end_loop:
  pop i
  ret

  .size matrix_multiply_SSE4_1_01, .-matrix_multiply_SSE4_1_01
