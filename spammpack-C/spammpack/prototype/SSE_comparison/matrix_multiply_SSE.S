# C API:
#
# void
# matrix_multiply_SSE (const unsigned int N, float *A, float *B, float *C);

#define N %rdi
#define A %rsi
#define B %rdx
#define C %rcx

#define i %rax

  .text
  .align 256
  .global matrix_multiply_SSE
  .type matrix_multiply_SSE, @function

matrix_multiply_SSE:

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

  # Calculate C(1,:).
  movaps 0x000(A), %xmm8
  movaps 0x010(A), %xmm9
  movaps 0x020(A), %xmm10
  mulps %xmm4, %xmm8
  mulps %xmm5, %xmm9
  addps %xmm8, %xmm0
  movaps 0x030(A), %xmm11
  mulps %xmm6, %xmm10
  addps %xmm9, %xmm0
  movaps 0x040(A), %xmm12
  mulps %xmm7, %xmm11
  addps %xmm10, %xmm0
  movaps 0x050(A), %xmm13
  mulps %xmm4, %xmm12
  addps %xmm11, %xmm0
  movaps 0x060(A), %xmm14
  mulps %xmm5, %xmm13
  addps %xmm12, %xmm1
  movaps 0x070(A), %xmm15
  mulps %xmm6, %xmm14
  addps %xmm13, %xmm1
  movaps 0x080(A), %xmm8
  mulps %xmm7, %xmm15
  addps %xmm14, %xmm1
  movaps 0x090(A), %xmm9
  mulps %xmm4, %xmm8
  addps %xmm15, %xmm1
  movaps 0x0a0(A), %xmm10
  mulps %xmm5, %xmm9
  addps %xmm8, %xmm2
  movaps 0x0b0(A), %xmm11
  mulps %xmm6, %xmm10
  addps %xmm9, %xmm2
  movaps 0x0c0(A), %xmm12
  mulps %xmm7, %xmm11
  addps %xmm10, %xmm2
  movaps 0x0d0(A), %xmm13
  mulps %xmm4, %xmm12
  addps %xmm11, %xmm2
  movaps 0x0e0(A), %xmm14
  mulps %xmm5, %xmm13
  addps %xmm12, %xmm3
  movaps 0x0f0(A), %xmm15
  mulps %xmm6, %xmm14
  addps %xmm13, %xmm3
  mulps %xmm7, %xmm15
  addps %xmm14, %xmm3
  addps %xmm15, %xmm3

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

  .size matrix_multiply_SSE, .-matrix_multiply_SSE
