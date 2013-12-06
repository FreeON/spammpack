#include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void
timer_start (struct timespec *const start_time)
{
  clock_gettime(CLOCKTYPE, start_time);
}

void
timer_stop (struct timespec *const end_time)
{
  clock_gettime(CLOCKTYPE, end_time);
}

double
timer_get (const struct timespec *const start_time,
    const struct timespec *const end_time)
{
  return end_time->tv_sec + end_time->tv_nsec/1.0e9
    - (start_time->tv_sec + start_time->tv_nsec/1.0e9);
}

void
work (void)
{
#pragma omp parallel for default(none)
  for(size_t i = 0; i < 100000000; i++)
  {
    double x = sin(i/1.0e-2);
    if(x < -2)
    {
      exit(1);
    }
  }
}

int
main ()
{
  struct timespec start_time;
  struct timespec end_time;

  double elapsed;

#ifdef _OPENMP
#pragma omp parallel
  printf("starting thread %d\n", omp_get_thread_num());
#endif

#ifdef _OPENMP
  printf("%d threads: ", omp_get_max_threads());
#endif
  timer_start(&start_time); work(); timer_stop(&end_time);
  elapsed = timer_get(&start_time, &end_time);
  printf("%e seconds", elapsed);

#ifdef _OPENMP
  printf(", %e seconds serial\n", elapsed*omp_get_max_threads());
#else
  printf(" in serial\n");
#endif

#ifdef _OPENMP
  omp_set_num_threads(1);
  printf("%d threads: ", omp_get_max_threads());
  timer_start(&start_time); work(); timer_stop(&end_time);
  elapsed = timer_get(&start_time, &end_time);
  printf("%e seconds\n", elapsed);
#endif
}
