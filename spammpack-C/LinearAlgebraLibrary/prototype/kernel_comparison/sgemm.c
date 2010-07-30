#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifdef HAVE_PAPI
#include <papi.h>
#endif

//#define EXTERNAL_BLAS
//#define NAIVE_4x4_KERNEL
//#define INLINE_SSE_4x4_KERNEL_1
#define EXTERNAL_SSE_4x4_KERNEL_1
//#define INLINE_SSE_4x4_KERNEL_2
//#define EXTERNAL_SSE_4x4_KERNEL_2

#if defined(EXTERNAL_SSE_4x4_KERNEL_1)
void
sgemm_kernel_1 (const unsigned int loops,
    const unsigned int N,
    float alpha, float *restrict A, float *restrict B,
    float beta, float *restrict C);
#endif

#if defined(EXTERNAL_SSE_4x4_KERNEL_2)
void
sgemm_kernel_2 (const unsigned int loops,
    const unsigned int N,
    float alpha, float *restrict A, float *restrict B,
    float beta, float *restrict C);
#endif

void
print_matrix (const unsigned int N, float *restrict A)
{
  unsigned int i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" %6.2f", A[i*N+j]);
    }
    printf("\n");
  }
}

void
multiply (const unsigned int loops,
    const unsigned int N,
    float alpha, float *restrict A, float *restrict B,
    float beta, float *restrict C)
{
  unsigned int loop;

#ifdef NAIVE_4x4_KERNEL
  unsigned int i, j, k;
#endif

#ifdef HAVE_PAPI
  int papi_events = PAPI_NULL;
  float papi_rtime;
  float papi_ptime;
  long_long papi_flpins;
  float papi_mflops;
  long long papi_ins;
  float papi_ipc;
#endif

  for (loop = 0; loop < loops; loop++)
  {
#ifdef HAVE_PAPI
    PAPI_create_eventset(&papi_events);
    PAPI_add_event(papi_events, PAPI_TOT_INS);
    PAPI_start(papi_events);
#endif

#if defined(EXTERNAL_BLAS)
    /* Use a standard external sgemm() from BLAS. */
    sgemm_("N", "N", &N, &N, &N, &alpha, A, &N, B, &N, &beta, C, &N);

#elif defined(NAIVE_4x4_KERNEL)
    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j)
      {
        C[i*4+j] *= beta;
      }
    }

    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j) {
        for (k = 0; k < 4; ++k)
        {
          C[i*4+j] += alpha*A[i*4+k]*B[k*4+j];
        }
      }
    }

#elif defined(INLINE_SSE_4x4_KERNEL_1)
    __asm__(
        /* Copy A into buffer. We will copy each matrix element of A 4 times
         * and arrange the elements in the buffer according to the following
         * storage layout:
         *
         * A_{11} A_{12} A_{13} A_{14}
         * A_{21} A_{22} A_{23} ...
         *
         * Copied to buffer as:
         *
         * A_{11} A_{11} A_{11} A_{11}
         * A_{12} A_{12} A_{12} A_{12}
         * A_{13} A_{13} ...
         */

        /* Make space for the buffer for B on the stack. */
        "movq %%rsp, %%r11\n\t" /* Save old stack. */
        "subq $(64 * 4 + 2 * 4), %%rsp\n\t"
        "andq $-4096, %%rsp\n\t" /* Align stack on a page boundary. */

        /* First we copy the first double word to the other double words in
         * xmm0 and xmm1. Unfortunately we have to go through memory, i.e. in
         * the asm inline we get alpha and beta in a normal register, but need
         * to move that into an SSE register. We do that by going to stack
         * first and then from stack to SSE. */
        "mov %3, 64*4(%%rsp)\n\t"
        "mov %4, 65*4(%%rsp)\n\t"

        "movss 64*4(%%rsp), %%xmm0\n\t"
        "movss 65*4(%%rsp), %%xmm1\n\t"

        "pshufd $0x00, %%xmm0, %%xmm0\n\t"
        "pshufd $0x00, %%xmm1, %%xmm1\n\t"

        /* Apply beta to C. */
        "movaps 0*4*4(%2), %%xmm2\n\t"
        "movaps 1*4*4(%2), %%xmm3\n\t"
        "movaps 2*4*4(%2), %%xmm4\n\t"
        "movaps 3*4*4(%2), %%xmm5\n\t"

        "mulps %%xmm1, %%xmm2\n\t"
        "mulps %%xmm1, %%xmm3\n\t"
        "mulps %%xmm1, %%xmm4\n\t"
        "mulps %%xmm1, %%xmm5\n\t"

        "movaps %%xmm2, 0*4*4(%2)\n\t"
        "movaps %%xmm3, 1*4*4(%2)\n\t"
        "movaps %%xmm4, 2*4*4(%2)\n\t"
        "movaps %%xmm5, 3*4*4(%2)\n\t"

        /* Copy matrix elements of A into SSE registers. */
        "movaps 0*4*4(%0), %%xmm3\n\t"
        "movaps 1*4*4(%0), %%xmm7\n\t"
        "movaps 2*4*4(%0), %%xmm11\n\t"
        "movaps 3*4*4(%0), %%xmm15\n\t"

        /* Apply alpha to A. */
        "mulps %%xmm0, %%xmm3\n\t"
        "mulps %%xmm0, %%xmm7\n\t"
        "mulps %%xmm0, %%xmm11\n\t"
        "mulps %%xmm0, %%xmm15\n\t"

        /* Dilate A. */
        "pshufd $0x00, %%xmm3, %%xmm0\n\t"
        "pshufd $0x55, %%xmm3, %%xmm1\n\t"
        "pshufd $0xaa, %%xmm3, %%xmm2\n\t"
        "pshufd $0xff, %%xmm3, %%xmm3\n\t"

        "pshufd $0x00, %%xmm7, %%xmm4\n\t"
        "pshufd $0x55, %%xmm7, %%xmm5\n\t"
        "pshufd $0xaa, %%xmm7, %%xmm6\n\t"
        "pshufd $0xff, %%xmm7, %%xmm7\n\t"

        "pshufd $0x00, %%xmm11, %%xmm8\n\t"
        "pshufd $0x55, %%xmm11, %%xmm9\n\t"
        "pshufd $0xaa, %%xmm11, %%xmm10\n\t"
        "pshufd $0xff, %%xmm11, %%xmm11\n\t"

        "pshufd $0x00, %%xmm15, %%xmm12\n\t"
        "pshufd $0x55, %%xmm15, %%xmm13\n\t"
        "pshufd $0xaa, %%xmm15, %%xmm14\n\t"
        "pshufd $0xff, %%xmm15, %%xmm15\n\t"

        /* Copy reordered A elements into buffer. */
        "movaps %%xmm0,   0*4*4(%%rsp)\n\t" /* A_{11} */
        "movaps %%xmm1,   1*4*4(%%rsp)\n\t" /* A_{12} */
        "movaps %%xmm2,   2*4*4(%%rsp)\n\t" /* A_{13} */
        "movaps %%xmm3,   3*4*4(%%rsp)\n\t" /* A_{14} */

        "movaps %%xmm4,   4*4*4(%%rsp)\n\t" /* A_{21} */
        "movaps %%xmm5,   5*4*4(%%rsp)\n\t" /* A_{22} */
        "movaps %%xmm6,   6*4*4(%%rsp)\n\t" /* A_{23} */
        "movaps %%xmm7,   7*4*4(%%rsp)\n\t" /* A_{24} */

        "movaps %%xmm8,   8*4*4(%%rsp)\n\t" /* A_{31} */
        "movaps %%xmm9,   9*4*4(%%rsp)\n\t" /* A_{32} */
        "movaps %%xmm10, 10*4*4(%%rsp)\n\t" /* A_{33} */
        "movaps %%xmm11, 11*4*4(%%rsp)\n\t" /* A_{34} */

        "movaps %%xmm12, 12*4*4(%%rsp)\n\t" /* A_{41} */
        "movaps %%xmm13, 13*4*4(%%rsp)\n\t" /* A_{42} */
        "movaps %%xmm14, 14*4*4(%%rsp)\n\t" /* A_{43} */
        "movaps %%xmm15, 15*4*4(%%rsp)\n\t" /* A_{44} */

        /* Zero out result matrix. */
        "pxor %%xmm8,  %%xmm8\n\t"
        "pxor %%xmm9,  %%xmm9\n\t"
        "pxor %%xmm10, %%xmm10\n\t"
        "pxor %%xmm11, %%xmm11\n\t"

        /* Load elements. */
        "movaps  0*4*4(%1),    %%xmm0\n\t" /* B_{11} B_{12} B_{13} B_{14} */
        "movaps  0*4*4(%%rsp), %%xmm2\n\t" /* A_{11} A_{11} A_{11} A_{11} */
        "movaps  4*4*4(%%rsp), %%xmm3\n\t" /* A_{21} A_{21} A_{21} A_{21} */
        "movaps  8*4*4(%%rsp), %%xmm4\n\t" /* A_{31} A_{31} A_{31} A_{31} */
        "movaps 12*4*4(%%rsp), %%xmm5\n\t" /* A_{41} A_{41} A_{41} A_{41} */

        /* Multiply A*B. */
        "mulps %%xmm0, %%xmm2\n\t"
        "mulps %%xmm0, %%xmm3\n\t"
        "mulps %%xmm0, %%xmm4\n\t"
        "mulps %%xmm0, %%xmm5\n\t"

        /* And store in C. */
        "addps %%xmm2, %%xmm8\n\t"
        "addps %%xmm3, %%xmm9\n\t"
        "addps %%xmm4, %%xmm10\n\t"
        "addps %%xmm5, %%xmm11\n\t"

        /* Load elements. */
        "movaps  1*4*4(%1),    %%xmm0\n\t" /* B_{21} B_{22} B_{23} B_{24} */
        "movaps  1*4*4(%%rsp), %%xmm2\n\t" /* A_{12} A_{12} A_{12} A_{12} */
        "movaps  5*4*4(%%rsp), %%xmm3\n\t" /* A_{22} A_{22} A_{22} A_{22} */
        "movaps  9*4*4(%%rsp), %%xmm4\n\t" /* A_{32} A_{32} A_{32} A_{32} */
        "movaps 13*4*4(%%rsp), %%xmm5\n\t" /* A_{42} A_{42} A_{42} A_{42} */

        /* Multiply A*B. */
        "mulps %%xmm0, %%xmm2\n\t"
        "mulps %%xmm0, %%xmm3\n\t"
        "mulps %%xmm0, %%xmm4\n\t"
        "mulps %%xmm0, %%xmm5\n\t"

        /* And store in C. */
        "addps %%xmm2, %%xmm8\n\t"
        "addps %%xmm3, %%xmm9\n\t"
        "addps %%xmm4, %%xmm10\n\t"
        "addps %%xmm5, %%xmm11\n\t"

        /* Load elements. */
        "movaps  2*4*4(%1),    %%xmm0\n\t" /* B_{31} B_{32} B_{33} B_{34} */
        "movaps  2*4*4(%%rsp), %%xmm2\n\t" /* A_{13} A_{13} A_{13} A_{13} */
        "movaps  6*4*4(%%rsp), %%xmm3\n\t" /* A_{23} A_{23} A_{23} A_{23} */
        "movaps 10*4*4(%%rsp), %%xmm4\n\t" /* A_{33} A_{33} A_{33} A_{33} */
        "movaps 14*4*4(%%rsp), %%xmm5\n\t" /* A_{43} A_{43} A_{43} A_{43} */

        /* Multiply A*B. */
        "mulps %%xmm0, %%xmm2\n\t"
        "mulps %%xmm0, %%xmm3\n\t"
        "mulps %%xmm0, %%xmm4\n\t"
        "mulps %%xmm0, %%xmm5\n\t"

        /* And store in C. */
        "addps %%xmm2, %%xmm8\n\t"
        "addps %%xmm3, %%xmm9\n\t"
        "addps %%xmm4, %%xmm10\n\t"
        "addps %%xmm5, %%xmm11\n\t"

        /* Load elements. */
        "movaps  3*4*4(%1),    %%xmm0\n\t" /* B_{41} B_{42} B_{43} B_{44} */
        "movaps  3*4*4(%%rsp), %%xmm2\n\t" /* A_{14} A_{14} A_{14} A_{14} */
        "movaps  7*4*4(%%rsp), %%xmm3\n\t" /* A_{24} A_{24} A_{24} A_{24} */
        "movaps 11*4*4(%%rsp), %%xmm4\n\t" /* A_{34} A_{34} A_{34} A_{34} */
        "movaps 15*4*4(%%rsp), %%xmm5\n\t" /* A_{44} A_{44} A_{44} A_{44} */

        /* Multiply A*B. */
        "mulps %%xmm0, %%xmm2\n\t"
        "mulps %%xmm0, %%xmm3\n\t"
        "mulps %%xmm0, %%xmm4\n\t"
        "mulps %%xmm0, %%xmm5\n\t"

        /* And store in C. */
        "addps %%xmm2, %%xmm8\n\t"
        "addps %%xmm3, %%xmm9\n\t"
        "addps %%xmm4, %%xmm10\n\t"
        "addps %%xmm5, %%xmm11\n\t"

        /* Add C from registers back into memory. */
        "movaps 0*4*4(%2), %%xmm0\n\t"
        "movaps 1*4*4(%2), %%xmm1\n\t"
        "movaps 2*4*4(%2), %%xmm2\n\t"
        "movaps 3*4*4(%2), %%xmm3\n\t"

        "addps %%xmm0, %%xmm8\n\t"
        "addps %%xmm1, %%xmm9\n\t"
        "addps %%xmm2, %%xmm10\n\t"
        "addps %%xmm3, %%xmm11\n\t"

        "movaps %%xmm8,  0*4*4(%2)\n\t"
        "movaps %%xmm9,  1*4*4(%2)\n\t"
        "movaps %%xmm10, 2*4*4(%2)\n\t"
        "movaps %%xmm11, 3*4*4(%2)\n\t"

        /* Restore stack. */
        "movq %%r11, %%rsp\n\t"

        : /* The output. */
        : /* The input. */
        "r" (A), "r" (B), "r" (C), "r" (alpha), "r" (beta)
        : /* Clobbered registers. */
        "r11", "xmm0", "xmm1", "xmm2", "xmm3", "xmm4", "xmm5", "xmm6", "xmm7",
        "xmm8", "xmm9", "xmm10", "xmm11", "xmm12", "xmm13", "xmm14", "xmm15"
           );

#elif defined(INLINE_SSE_4x4_KERNEL_2)
    __asm__(
        ".intel_syntax noprefix\n\t"

        "movaps   xmm6, [%0]\n\t"             /* Load A(1,1:4) into xmm0. */

        "movaps   xmm0, [%1]\n\t"             /* Load B(1,1:4) into xmm1. */
        "mulps    xmm0, xmm6\n\t"             /* Multiply A(1,1:4).*B(1,1:4) and store result in xmm1. */

        "movaps   xmm1, [%1+1*4*4]\n\t"       /* Load B(2,1:4) into xmm3. */
        "mulps    xmm1, xmm6\n\t"             /* Multiply A(1,1:4).*B(2,1:4) and store result in xmm3. */

        "movaps   xmm2, [%1+2*4*4]\n\t"       /* Load B(3,1:4) into xmm5. */
        "mulps    xmm2, xmm6\n\t"             /* Multiply A(1,1:4).*B(3,1:4) and store result in xmm5. */

        "movaps   xmm3, [%1+3*4*4]\n\t"       /* Load B(4,1:4) into xmm7. */
        "mulps    xmm3, xmm6\n\t"             /* Multiply A(1,1:4).*B(4,1:4) and store result in xmm7. */

        "haddps   xmm0, xmm1\n\t"
        "haddps   xmm2, xmm3\n\t"
        "haddps   xmm0, xmm2\n\t"
        "movaps   [%2], xmm0\n\t"

        "movaps   xmm6, [%0+1*4*4]\n\t"       /* Load A(2,1:4) into xmm0. */

        "movaps   xmm0, [%1]\n\t"             /* Load B(1,1:4) into xmm1. */
        "mulps    xmm0, xmm6\n\t"             /* Multiply A(1,1:4).*B(1,1:4) and store result in xmm1. */

        "movaps   xmm1, [%1+1*4*4]\n\t"       /* Load B(2,1:4) into xmm3. */
        "mulps    xmm1, xmm6\n\t"             /* Multiply A(1,1:4).*B(2,1:4) and store result in xmm3. */

        "movaps   xmm2, [%1+2*4*4]\n\t"       /* Load B(3,1:4) into xmm5. */
        "mulps    xmm2, xmm6\n\t"             /* Multiply A(1,1:4).*B(3,1:4) and store result in xmm5. */

        "movaps   xmm3, [%1+3*4*4]\n\t"       /* Load B(4,1:4) into xmm7. */
        "mulps    xmm3, xmm6\n\t"             /* Multiply A(1,1:4).*B(4,1:4) and store result in xmm7. */

        "haddps   xmm0, xmm1\n\t"
        "haddps   xmm2, xmm3\n\t"
        "haddps   xmm0, xmm2\n\t"
        "movaps   [%2+4*4], xmm0\n\t"

        "movaps   xmm6, [%0+2*4*4]\n\t"       /* Load A(3,1:4) into xmm0. */

        "movaps   xmm0, [%1]\n\t"             /* Load B(1,1:4) into xmm1. */
        "mulps    xmm0, xmm6\n\t"             /* Multiply A(1,1:4).*B(1,1:4) and store result in xmm1. */

        "movaps   xmm1, [%1+1*4*4]\n\t"       /* Load B(2,1:4) into xmm3. */
        "mulps    xmm1, xmm6\n\t"             /* Multiply A(1,1:4).*B(2,1:4) and store result in xmm3. */

        "movaps   xmm2, [%1+2*4*4]\n\t"       /* Load B(3,1:4) into xmm5. */
        "mulps    xmm2, xmm6\n\t"             /* Multiply A(1,1:4).*B(3,1:4) and store result in xmm5. */

        "movaps   xmm3, [%1+3*4*4]\n\t"       /* Load B(4,1:4) into xmm7. */
        "mulps    xmm3, xmm6\n\t"             /* Multiply A(1,1:4).*B(4,1:4) and store result in xmm7. */

        "haddps   xmm0, xmm1\n\t"
        "haddps   xmm2, xmm3\n\t"
        "haddps   xmm0, xmm2\n\t"
        "movaps   [%2+2*4*4], xmm0\n\t"

        "movaps   xmm6, [%0+3*4*4]\n\t"       /* Load A(4,1:4) into xmm0. */

        "movaps   xmm0, [%1]\n\t"             /* Load B(1,1:4) into xmm1. */
        "mulps    xmm0, xmm6\n\t"             /* Multiply A(1,1:4).*B(1,1:4) and store result in xmm1. */

        "movaps   xmm1, [%1+1*4*4]\n\t"       /* Load B(2,1:4) into xmm3. */
        "mulps    xmm1, xmm6\n\t"             /* Multiply A(1,1:4).*B(2,1:4) and store result in xmm3. */

        "movaps   xmm2, [%1+2*4*4]\n\t"       /* Load B(3,1:4) into xmm5. */
        "mulps    xmm2, xmm6\n\t"             /* Multiply A(1,1:4).*B(3,1:4) and store result in xmm5. */

        "movaps   xmm3, [%1+3*4*4]\n\t"       /* Load B(4,1:4) into xmm7. */
        "mulps    xmm3, xmm6\n\t"             /* Multiply A(1,1:4).*B(4,1:4) and store result in xmm7. */

        "haddps   xmm0, xmm1\n\t"
        "haddps   xmm2, xmm3\n\t"
        "haddps   xmm0, xmm2\n\t"
        "movaps   [%2+3*4*4], xmm0\n\t"

        ".att_syntax prefix\n\t"
        : /* The output. */
        : /* The input. */
        "r" (A), "r" (B), "r" (C), "r" (alpha), "r" (beta)
        : /* Clobbered registers. */
          "xmm0", "xmm1", "xmm2", "xmm3", "xmm4",
          "xmm5", "xmm6"
           );
#elif defined(EXTERNAL_SSE_4x4_KERNEL_1)
    sgemm_kernel_1(loops, N, alpha, A, B, beta, C);
#elif defined(EXTERNAL_SSE_4x4_KERNEL_2)
    sgemm_kernel_2(loops, N, alpha, A, B, beta, C);
#endif

#ifdef HAVE_PAPI
    PAPI_stop(papi_events, &papi_ins);
    printf("PAPI (in loop) - %lli total instructions\n", papi_ins);
    PAPI_cleanup_eventset(papi_events);
    PAPI_destroy_eventset(&papi_events);
#endif
  }
}

int
main (int argc, char **argv)
{
  int loops = 1;
  unsigned int N = 1024;

  int i, j, k, loop;

  float max_diff = 0;

  float alpha = 1.2;
  float beta = 0.5;
  float *A, *B, *C, *D;

  struct timeval start, stop;
  struct rusage rusage_start, rusage_stop;
  double walltime, usertime, systime, flops;

  int verify = 0;

  int parse;
  int longindex;
  char *short_options = "hN:l:v";
  struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "N", required_argument, NULL, 'N' },
    { "loops", required_argument, NULL, 'l' },
    { "verify", no_argument, NULL, 'v' },
    { NULL, 0, NULL, 0 }
  };

  /* Read command line. */
  while ((parse = getopt_long(argc, argv, short_options, long_options, &longindex)) != -1)
  {
    switch (parse)
    {
      case 'h':
        printf("Usage:\n");
        printf("\n");
        printf("-h            This help\n");
        printf("-N N          Use NxN matrix blocks\n");
        printf("--loops N     Repeat each multiply N times\n");
        printf("--verify      Verify result\n");
        return 0;
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'l':
        loops = strtol(optarg, NULL, 10);
        break;

      case 'v':
        verify = 1;
        break;

      default:
        printf("unknown command line argument\n");
        return -1;
        break;
    }
  }

#if defined(EXTERNAL_BLAS)
  printf("external blas\n");
#elif defined(NAIVE_4x4_KERNEL)
  printf("naive 4x4 kernel\n");
#elif defined(INLINE_SSE_4x4_KERNEL_1)
  printf("inline sgemm_kernel_1\n");
#elif defined(EXTERNAL_SSE_4x4_KERNEL_1)
  printf("external sgemm_kernel_1\n");
#elif defined(INLINE_SSE_4x4_KERNEL_2)
  printf("inline sgemm_kernel_2\n");
#elif defined(EXTERNAL_SSE_4x4_KERNEL_2)
  printf("external sgemm_kernel_2\n");
#endif

  printf("alpha = %f, beta = %f\n", alpha, beta);

  if (verify && loops > 1)
  {
    printf("can not verify for loops > 1\n");
    verify = 0;
  }

#if defined(NAIVE_4x4_KERNEL) || defined(INLINE_SSE_4x4_KERNEL_1) || defined(INLINE_SSE_4x4_KERNEL_2) || defined(EXTERNAL_SSE_4x4_KERNEL_1) || defined(EXTERNAL_SSE_4x4_KERNEL_2)
  if (N != 4)
  {
    printf("This kernel only works with 4x4 matrices\n");
    N = 4;
  }
#endif

#ifdef HAVE_POSIX_MEMALIGN
  if (posix_memalign((void**) &A, 64, sizeof(float)*N*N) != 0)
  {
    printf("error allocating matrix memory for A\n");
    exit(1);
  }

  if (posix_memalign((void**) &B, 64, sizeof(float)*N*N) != 0)
  {
    printf("error allocating matrix memory for B\n");
    exit(1);
  }

  if (posix_memalign((void**) &C, 64, sizeof(float)*N*N) != 0)
  {
    printf("error allocating matrix memory for C\n");
    exit(1);
  }

  if (posix_memalign((void**) &D, 64, sizeof(float)*N*N) != 0)
  {
    printf("error allocating matrix memory for D\n");
    exit(1);
  }
#else
  A = (float*) malloc(sizeof(float)*N*N);
  B = (float*) malloc(sizeof(float)*N*N);
  C = (float*) malloc(sizeof(float)*N*N);
  D = (float*) malloc(sizeof(float)*N*N);
#endif

#ifdef HAVE_PAPI
  /* Do some PAPI. */
  if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
  {
    printf("can not initialize PAPI\n");
    exit(1);
  }
#endif

  printf("filling matrices with random matrix elements\n");
  for (i = 0; i < N; ++i)
    for (j = 0; j < N; ++j) {
    {
      //A[i*N+j] = rand()/(double) RAND_MAX;
      //B[i*N+j] = rand()/(double) RAND_MAX;
      //C[i*N+j] = rand()/(double) RAND_MAX;

      //A[i*N+j] = i*N+j+1;
      //B[i*N+j] = i*N+j+1+N*N;
      //C[i*N+j] = i*N+j+1+2*N*N;

      A[i*N+j] = i*N+j+1;
      B[i*N+j] = i*N+j+1;
      C[i*N+j] = i*N+j+1;
      D[i*N+j] = C[i*N+j];
    }
  }

  printf("looping %u times\n", loops);
  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);
  multiply(loops, N, alpha, A, B, beta, C);
  getrusage(RUSAGE_SELF, &rusage_stop);
  gettimeofday(&stop, NULL);
  walltime = (stop.tv_sec-start.tv_sec+(stop.tv_usec-start.tv_usec)/1.0e6)/loops;
  usertime = ((rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6)/loops;
  systime = ((rusage_stop.ru_stime.tv_sec-rusage_start.ru_stime.tv_sec)+(rusage_stop.ru_stime.tv_usec-rusage_start.ru_stime.tv_usec)/(double) 1e6)/loops;
  flops = ((double) N)*((double) N)*(2.*N+1.)/usertime;
  printf("%ix%i matrix (repeated %i times)\n", N, N, loops);
  if (flops < 1000*1000*1000)
  {
    printf("performance: total walltime %f s, usertime %f s, systime %f s, per iteration walltime %e s, usertime %e s, systime %e s, %1.2f Mflop/s\n",
        walltime*loops, usertime*loops, systime*loops,
        walltime, usertime, systime,
        flops/1000./1000.);
  }

  else
  {
    printf("performance: total walltime %f s, usertime %f s, systime %f s, per iteration walltime %e s, usertime %e s, systime %e s, %1.2f Gflop/s\n",
        walltime*loops, usertime*loops, systime*loops,
        walltime, usertime, systime,
        flops/1000./1000./1000.);
  }

  /* Check result. */
  if (verify)
  {
    printf("verifying result\n");
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++)
      {
        D[i*N+j] *= beta;
      }
    }

    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        for (k = 0; k < N; k++)
        {
          D[i*N+j] += alpha*A[i*N+k]*B[k*N+j];
        }
      }
    }

    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++)
      {
        if (fabs(C[i*N+j]-D[i*N+j]) > max_diff)
        {
          max_diff = fabs(C[i*N+j]-D[i*N+j]);
        }
      }
    }

    if (max_diff > 1e-14)
    {
      printf("C =\n");
      print_matrix(N, C);

      printf("D =\n");
      print_matrix(N, D);

      printf("max diff = %e\n", max_diff);
    }

    else
    {
      printf("result correct\n");
    }
  }

  free(A);
  free(B);
  free(C);
}
