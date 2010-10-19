#define loops     %rdi
#define N         %rsi
#define alpha     %xmm0
#define A         %rdx
#define B         %rcx
#define beta      %xmm1
#define C         %r8
#define OLD_STACK %r11

  .text
  .align 256
  .global sgemm_kernel_2
  .type sgemm_kernel_2, @function

sgemm_kernel_2:

  .intel_syntax noprefix

  movaps   xmm6, [A]              /* Load A(1,1:4) into xmm0. */

  movaps   xmm0, [B]              /* Load B(1,1:4) into xmm1. */
  mulps    xmm0, xmm6             /* Multiply A(1,1:4).*B(1,1:4) and store result in xmm1. */

  movaps   xmm1, [B+1*4*4]        /* Load B(2,1:4) into xmm3. */
  mulps    xmm1, xmm6             /* Multiply A(1,1:4).*B(2,1:4) and store result in xmm3. */

  movaps   xmm2, [B+2*4*4]        /* Load B(3,1:4) into xmm5. */
  mulps    xmm2, xmm6             /* Multiply A(1,1:4).*B(3,1:4) and store result in xmm5. */

  movaps   xmm3, [B+3*4*4]        /* Load B(4,1:4) into xmm7. */
  mulps    xmm3, xmm6             /* Multiply A(1,1:4).*B(4,1:4) and store result in xmm7. */

  haddps   xmm0, xmm1
  haddps   xmm2, xmm3
  haddps   xmm0, xmm2
  movaps   [C], xmm0

  movaps   xmm6, [A+1*4*4]        /* Load A(2,1:4) into xmm0. */

  movaps   xmm0, [B]              /* Load B(1,1:4) into xmm1. */
  mulps    xmm0, xmm6             /* Multiply A(1,1:4).*B(1,1:4) and store result in xmm1. */

  movaps   xmm1, [B+1*4*4]        /* Load B(2,1:4) into xmm3. */
  mulps    xmm1, xmm6             /* Multiply A(1,1:4).*B(2,1:4) and store result in xmm3. */

  movaps   xmm2, [B+2*4*4]        /* Load B(3,1:4) into xmm5. */
  mulps    xmm2, xmm6             /* Multiply A(1,1:4).*B(3,1:4) and store result in xmm5. */

  movaps   xmm3, [B+3*4*4]        /* Load B(4,1:4) into xmm7. */
  mulps    xmm3, xmm6             /* Multiply A(1,1:4).*B(4,1:4) and store result in xmm7. */

  haddps   xmm0, xmm1
  haddps   xmm2, xmm3
  haddps   xmm0, xmm2
  movaps   [C+4*4], xmm0

  movaps   xmm6, [A+2*4*4]        /* Load A(3,1:4) into xmm0. */

  movaps   xmm0, [B]              /* Load B(1,1:4) into xmm1. */
  mulps    xmm0, xmm6             /* Multiply A(1,1:4).*B(1,1:4) and store result in xmm1. */

  movaps   xmm1, [B+1*4*4]        /* Load B(2,1:4) into xmm3. */
  mulps    xmm1, xmm6             /* Multiply A(1,1:4).*B(2,1:4) and store result in xmm3. */

  movaps   xmm2, [B+2*4*4]        /* Load B(3,1:4) into xmm5. */
  mulps    xmm2, xmm6             /* Multiply A(1,1:4).*B(3,1:4) and store result in xmm5. */

  movaps   xmm3, [B+3*4*4]        /* Load B(4,1:4) into xmm7. */
  mulps    xmm3, xmm6             /* Multiply A(1,1:4).*B(4,1:4) and store result in xmm7. */

  haddps   xmm0, xmm1
  haddps   xmm2, xmm3
  haddps   xmm0, xmm2
  movaps   [C+2*4*4], xmm0

  movaps   xmm6, [A+3*4*4]        /* Load A(4,1:4) into xmm0. */

  movaps   xmm0, [B]              /* Load B(1,1:4) into xmm1. */
  mulps    xmm0, xmm6             /* Multiply A(1,1:4).*B(1,1:4) and store result in xmm1. */

  movaps   xmm1, [B+1*4*4]        /* Load B(2,1:4) into xmm3. */
  mulps    xmm1, xmm6             /* Multiply A(1,1:4).*B(2,1:4) and store result in xmm3. */

  movaps   xmm2, [B+2*4*4]        /* Load B(3,1:4) into xmm5. */
  mulps    xmm2, xmm6             /* Multiply A(1,1:4).*B(3,1:4) and store result in xmm5. */

  movaps   xmm3, [B+3*4*4]        /* Load B(4,1:4) into xmm7. */
  mulps    xmm3, xmm6             /* Multiply A(1,1:4).*B(4,1:4) and store result in xmm7. */

  haddps   xmm0, xmm1
  haddps   xmm2, xmm3
  haddps   xmm0, xmm2
  movaps   [C+3*4*4], xmm0

  .att_syntax prefix

  /* Done. */
  ret

  .size sgemm_kernel_2, .-sgemm_kernel_2
