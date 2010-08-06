#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifdef HAVE_PAPI
#include <papi.h>
#endif

#define CACHELINE_SIZE 64

int
main (int argc, char **argv)
{
  unsigned long long loop;
  unsigned long long loops = 1;
  unsigned long long N = 1024*1024*1024;
  unsigned long long d = 1;

  struct timeval start, stop;
  struct rusage rusage_start, rusage_stop;
  double walltime, usertime, systime;

#ifdef HAVE_PAPI
  int papi_events = PAPI_NULL;
  long long *papi_values;
#endif

  float *A;

  int parse;
  int longindex;
  char *short_options = "hN:d:l:";
  struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "N", required_argument, NULL, 'N' },
    { "d", required_argument, NULL, 'd' },
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
        printf("-h            This help\n");
        printf("-N N          Use N byte buffer\n");
        printf("-d d          access 2 bytes d bytes apart\n");
        printf("--loops N     Repeat each access test N times\n");
        return 0;
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'd':
        d = strtol(optarg, NULL, 10);
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

  if (2*d > N)
  {
    printf("distance has to be at most N/2 (%llu)\n", N/2);
    exit(1);
  }

#ifdef HAVE_POSIX_MEMALIGN
  if (posix_memalign((void**) &A, CACHELINE_SIZE, sizeof(float)*N) != 0)
  {
    printf("error allocating A\n");
    exit(1);
  }
#else
  if ((A = (struct float*) malloc(sizeof(float)*N)) == NULL)
  {
    printf("error allocating A\n");
    exit(1);
  }
#endif

#ifdef HAVE_PAPI
  /* Do some PAPI. */
  if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
  {
    printf("can not initialize PAPI\n");
    exit(1);
  }
#endif

  if (N < 1024) printf("buffer is %llu bytes\n", N*sizeof(float));
  else if (N < 1024*1024) printf("buffer is %1.1f kB\n", N*sizeof(float)/1024.);
  else if (N < 1024*1024*1024) printf("buffer is %1.1f MB\n", N*sizeof(float)/1024./1024.);
  else printf("buffer is %1.1f GB\n", N*sizeof(float)/1024./1024./1024.);
  printf("accessing A[0] = A[%llu]*A[%llu] --> %p = %p * %p\n", d, 2*d, (void*) &A[0], (void*) &A[d], (void*) &A[2*d]);

  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);

#ifdef HAVE_PAPI
  papi_values = (long long*) malloc(sizeof(long long)*6);
  PAPI_set_granularity(PAPI_GRN_MIN);
  PAPI_create_eventset(&papi_events);
  PAPI_add_event(papi_events, PAPI_TOT_INS);
  PAPI_add_event(papi_events, PAPI_TOT_CYC);
  PAPI_add_event(papi_events, PAPI_RES_STL);
  PAPI_add_event(papi_events, PAPI_L1_TCM);
  PAPI_add_event(papi_events, PAPI_L2_TCM);
  PAPI_add_event(papi_events, PAPI_TLB_DM);
  PAPI_start(papi_events);
#endif

  for (loop = 0; loop < loops; loop++)
  {
    A[0] = A[d]*A[2*d];
  }

#ifdef HAVE_PAPI
  PAPI_stop(papi_events, papi_values);
  printf("[PAPI] %lli total instructions\n",    papi_values[0]);
  printf("[PAPI] %lli total cycles\n",          papi_values[1]);
  printf("[PAPI] %lli cycles stalled\n",        papi_values[2]);
  printf("[PAPI] %lli total L1 misses\n",       papi_values[3]);
  printf("[PAPI] %lli total L2 misses\n",       papi_values[4]);
  printf("[PAPI] %lli total data TLB misses\n", papi_values[5]);
  PAPI_cleanup_eventset(papi_events);
  PAPI_destroy_eventset(&papi_events);
  free(papi_values);
#endif

  getrusage(RUSAGE_SELF, &rusage_stop);
  gettimeofday(&stop, NULL);
  walltime = (stop.tv_sec-start.tv_sec+(stop.tv_usec-start.tv_usec)/1.0e6)/loops;
  usertime = ((rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6)/loops;
  systime = ((rusage_stop.ru_stime.tv_sec-rusage_start.ru_stime.tv_sec)+(rusage_stop.ru_stime.tv_usec-rusage_start.ru_stime.tv_usec)/(double) 1e6)/loops;
  printf("performance: total walltime %f s, usertime %f s, systime %f s, per iteration walltime %e s, usertime %e s, systime %e s\n",
      walltime*loops, usertime*loops, systime*loops, walltime, usertime, systime);

  free(A);
}
