#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/resource.h>

#include <xmmintrin.h>

//#define EXTERNAL_BLAS
//#define NAIVE_4x4_KERNEL
#define SSE_4x4_KERNEL

void
multiply (const unsigned int loops,
    const unsigned int N,
    float alpha, float *A, float *B,
    float beta, float *C)
{
  unsigned int loop;
  unsigned int i, j, k;

  float *D = (float*) malloc(sizeof(float)*4*4);
  float max_diff = 0;

  for (loop = 0; loop < loops; loop++)
  {

    /* Use a standard external sgemm() from BLAS. */
#ifdef EXTERNAL_BLAS
    sgemm_("N", "N", &N, &N, &N, &alpha, A, &N, B, &N, &beta, C, &N);
#endif

#ifdef NAIVE_4x4_KERNEL
    for (i = 0; i < 4; ++i) {
      for (j = 0; j < 4; ++j) {
        for (k = 0; k < 4; ++k)
        {
          C[i*4+j] += A[i*4+k]*B[k*4+j];
        }
      }
    }
#endif

#ifdef SSE_4x4_KERNEL
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
#endif

  }
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

#if defined(NAIVE_4x4_KERNEL) || defined(SSE_4x4_KERNEL)
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

  for (j = 0; j < N; ++j) {
    for (i = 0; i < N; ++i)
    {
      A[i+j*N] = rand()/(double) RAND_MAX;
      B[i+j*N] = rand()/(double) RAND_MAX;
      C[i+j*N] = rand()/(double) RAND_MAX;
    }
  }

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
