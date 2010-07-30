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
//#define SSE_4x4_KERNEL
//#define BETTER_SSE_4x4_KERNEL
#define EXTERNAL_SSE_4x4_KERNEL

#if defined(EXTERNAL_SSE_4x4_KERNEL)
void
sgemm_kernel_1 (const unsigned int loops,
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
  unsigned int i, j, k;

  float *D = (float*) malloc(sizeof(float)*4*4);
  float max_diff = 0;

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
      for (j = 0; j < 4; ++j) {
        for (k = 0; k < 4; ++k)
        {
          C[i*4+j] += A[i*4+k]*B[k*4+j];
        }
      }
    }

#elif defined(SSE_4x4_KERNEL)
    //printf("A(1,1:4) = [ %f %f %f %f ]\n", A[0], A[1], A[2], A[3]);
    //printf("B(1,1:4) = [ %f %f %f %f ]\n", B[0], B[1], B[2], B[3]);

    /* Assume B is transposed. */
    //for (i = 0; i < 4; i++) {
    //  for (j = 0; j < 4; j++)
    //  {
    //    D[i*4+j] = 0;
    //    for (k = 0; k < 4; k++)
    //    {
    //      D[i*4+j] += A[i*4+k]*B[j*4+k];
    //    }
    //  }
    //}

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

    //printf("C =\n");
    //for (i = 0; i < 4; i++) {
    //  for (j = 0; j < 4; j++)
    //  {
    //    printf(" %f", C[i*4+j]);
    //  }
    //  printf("\n");
    //}

    //printf("D =\n");
    //for (i = 0; i < 4; i++) {
    //  for (j = 0; j < 4; j++)
    //  {
    //    printf(" %f", D[i*4+j]);
    //  }
    //  printf("\n");
    //}

    //for (i = 0; i < 4; i++) {
    //  for (j = 0; j < 4; j++)
    //  {
    //    if (fabs(C[i*4+j]-D[i*4+j]) > max_diff)
    //    {
    //      max_diff = fabs(C[i*4+j]-D[i*4+j]);
    //    }
    //  }
    //}
    //printf("max diff = %f\n", max_diff);

    //exit(0);

#elif defined(BETTER_SSE_4x4_KERNEL)
    __asm__(
        /* Copy B into buffer. We will copy each matrix element of B 4 times
         * and arrange the elements in the buffer according to the following
         * storage layout:
         *
         * B_{11} B_{11} B_{11} B_{11}
         * B_{12} B_{12} B_{12} B_{12}
         * B_{13} B_{13} ...
         */

        /* Make space for the buffer on the stack. */
        "subq $(64 * 4), %%rsp\n\t"

        /* Copy matrix elements into buffer. */
        "movaps 0*4*4(%1), %%xmm3\n\t"
        "movaps 1*4*4(%1), %%xmm7\n\t"
        "movaps 2*4*4(%1), %%xmm11\n\t"
        "movaps 3*4*4(%1), %%xmm15\n\t"

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

        "movaps %%xmm0,   0*4*4(%%rsp)\n\t"
        "movaps %%xmm1,   1*4*4(%%rsp)\n\t"
        "movaps %%xmm2,   2*4*4(%%rsp)\n\t"
        "movaps %%xmm3,   3*4*4(%%rsp)\n\t"

        "movaps %%xmm4,   4*4*4(%%rsp)\n\t"
        "movaps %%xmm5,   5*4*4(%%rsp)\n\t"
        "movaps %%xmm6,   6*4*4(%%rsp)\n\t"
        "movaps %%xmm7,   7*4*4(%%rsp)\n\t"

        "movaps %%xmm8,   8*4*4(%%rsp)\n\t"
        "movaps %%xmm9,   9*4*4(%%rsp)\n\t"
        "movaps %%xmm10, 10*4*4(%%rsp)\n\t"
        "movaps %%xmm11, 11*4*4(%%rsp)\n\t"

        "movaps %%xmm12, 12*4*4(%%rsp)\n\t"
        "movaps %%xmm13, 13*4*4(%%rsp)\n\t"
        "movaps %%xmm14, 14*4*4(%%rsp)\n\t"
        "movaps %%xmm15, 15*4*4(%%rsp)\n\t"

        /* Remove buffer from stack. */
        "addq $(64 * 4), %%rsp\n\t"
        : /* The output. */
        : /* The input. */
        "r" (A), "r" (B), "r" (C), "r" (alpha), "r" (beta)
        : /* Clobbered registers. */
          "xmm0", "xmm1", "xmm2", "xmm3", "xmm4",
          "xmm5", "xmm6", "xmm7", "xmm8", "xmm9",
          "xmm10", "xmm11", "xmm12", "xmm13", "xmm14",
          "xmm15"
           );
#elif defined(EXTERNAL_SSE_4x4_KERNEL)
    sgemm_kernel_1(loops, N, alpha, A, B, beta, C);
#endif

#ifdef HAVE_PAPI
    PAPI_stop(papi_events, &papi_ins);
    printf("PAPI (in loop) - %lli total instructions\n", papi_ins);
    PAPI_cleanup_eventset(papi_events);
    PAPI_destroy_eventset(&papi_events);
#endif
  }

  free(D);
}

int
main (int argc, char **argv)
{
  int loops = 1;
  unsigned int N = 1024;

  int i, j, loop;

  float alpha = 1.2;
  float beta = 0.5;
  float *A, *B, *C;

  struct timeval start, stop;
  struct rusage rusage_start, rusage_stop;
  double walltime, usertime, systime, flops;

  int parse;
  int longindex;
  char *short_options = "hN:l:";
  struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "N", required_argument, NULL, 'N' },
    { "loops", required_argument, NULL, 'l' },
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
        printf("-h           This help\n");
        printf("-N N         Use NxN matrix blocks\n");
        printf("--loops N    Repeat each multiply N times\n");
        return 0;
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'l':
        loops = strtol(optarg, NULL, 10);
        break;

      default:
        printf("unknown command line argument\n");
        return -1;
        break;
    }
  }

#if defined(NAIVE_4x4_KERNEL) || defined(SSE_4x4_KERNEL) || defined(BETTER_SSE_4x4_KERNEL) | defined(EXTERNAL_SSE_4x4_KERNEL)
  if (N != 4)
  {
    printf("This kernel only works with 4x4 matrices\n");
    N = 4;
  }
#endif

  if (posix_memalign((void**) &A, 64, sizeof(float)*N*N) != 0)
  {
    printf("error allocating matrix memory\n");
    exit(1);
  }

  if (posix_memalign((void**) &B, 64, sizeof(float)*N*N) != 0)
  {
    printf("error allocating matrix memory\n");
    exit(1);
  }

  if (posix_memalign((void**) &C, 64, sizeof(float)*N*N) != 0)
  {
    printf("error allocating matrix memory\n");
    exit(1);
  }

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
  flops = ((double) N)*((double) N)*(2.*N+1.)/walltime;
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

  free(A);
  free(B);
  free(C);
}
