#include <stdlib.h>
#include <stdio.h>
#include <emmintrin.h>
#include <sys/time.h>
#include <sys/resource.h>

#define N 2048

double res[N][N]  __attribute__ ((aligned (64)));
double mul1[N][N] __attribute__ ((aligned (64)));
double mul2[N][N] __attribute__ ((aligned (64)));

#define CLS 64
#define SM (CLS / sizeof (double))

int
main (void)
{
  int i, i2, j, j2, k, k2;
  double *restrict rres;
  double *restrict rmul1;
  double *restrict rmul2;

  struct timeval start, stop;
  struct rusage rusage_start, rusage_stop;
  double walltime, usertime, systime, flops;

  // ... Initialize mul1 and mul2
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      mul1[i][j] = rand()/(double) RAND_MAX;
      mul2[i][j] = rand()/(double) RAND_MAX;
    }
  }

  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);

  for (i = 0; i < N; i += SM)
    for (j = 0; j < N; j += SM)
      for (k = 0; k < N; k += SM)
        for (i2 = 0, rres = &res[i][j], rmul1 = &mul1[i][k]; i2 < SM; ++i2, rres += N, rmul1 += N)
        {
          _mm_prefetch (&rmul1[8], _MM_HINT_NTA);
          for (k2 = 0, rmul2 = &mul2[k][j]; k2 < SM; ++k2, rmul2 += N)
          {
            __m128d m1d = _mm_load_sd (&rmul1[k2]);
            m1d = _mm_unpacklo_pd (m1d, m1d);
            for (j2 = 0; j2 < SM; j2 += 2)
            {
              __m128d m2 = _mm_load_pd (&rmul2[j2]);
              __m128d r2 = _mm_load_pd (&rres[j2]);
              _mm_store_pd (&rres[j2], _mm_add_pd (_mm_mul_pd (m2, m1d), r2));
            }
          }
        }

  getrusage(RUSAGE_SELF, &rusage_stop);
  gettimeofday(&stop, NULL);
  walltime = (stop.tv_sec-start.tv_sec+(stop.tv_usec-start.tv_usec)/1.0e6);
  usertime = ((rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6);
  systime = ((rusage_stop.ru_stime.tv_sec-rusage_start.ru_stime.tv_sec)+(rusage_stop.ru_stime.tv_usec-rusage_start.ru_stime.tv_usec)/(double) 1e6);
  flops = ((double) N)*((double) N)*(2.*N+1.)/usertime;
  printf("performance: total walltime %f s, usertime %f s, systime %f s, ", walltime, usertime, systime);
  if (flops < 1000*1000*1000)
  {
    printf("%1.2f Mflop/s\n", flops/1000./1000.);
  }

  else
  {
    printf("%1.2f Gflop/s\n", flops/1000./1000./1000.);
  }

  // ... use res matrix

  return 0;
}
