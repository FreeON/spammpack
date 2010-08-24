#include "config.h"
#include <emmintrin.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifdef HAVE_PAPI
#include <papi.h>
#endif

#define N_BLOCK 4

#ifdef HAVE_PAPI
void
print_papi_error (const char *message, const int errorcode)
{
  char error[1001];

  PAPI_perror(errorcode, error, 1000);
  printf("PAPI error: %s - %s\n", message, error);
  exit(1);
}
#endif

int
main (int argc, char **argv)
{
  unsigned long long loops = 1;
  unsigned long long N = 1;
  unsigned long long alignment = 4096;

  int result;

  unsigned long long loop;
  unsigned long long index;
  unsigned int i, j;

  struct timeval start, stop;
  struct rusage rusage_start, rusage_stop;
  double walltime, usertime, systime, flops;

#ifdef HAVE_PAPI
  int papi_events = PAPI_NULL;
  long long *papi_values;
  int papi_value_index;
  int papi_errorcode;

  int load_TOT_INS = 0;
  int load_TOT_CYC = 0;
  int load_RES_STL = 0;
  int load_L1_ICM = 0;
  int load_L1_DCM = 0;
  int load_L2_ICM = 0;
  int load_L2_DCM = 0;
  int load_TLB_IM = 0;
  int load_TLB_DM = 0;
#endif

  float *A, *B, *C, *stream;
  float *restrict A_block, *restrict B_block, *restrict C_block;

  int parse;
  int longindex;
  char *short_options = "hN:a:d:l:123456789";
  struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "N", required_argument, NULL, 'N' },
    { "align", required_argument, NULL, 'a' },
    { "d", required_argument, NULL, 'd' },
    { "loops", required_argument, NULL, 'l' },
#ifdef HAVE_PAPI
    { "TOT_INS", no_argument, NULL, '1' },
    { "TOT_CYC", no_argument, NULL, '2' },
    { "RES_STL", no_argument, NULL, '3' },
    { "L1_ICM", no_argument, NULL, '4' },
    { "L1_DCM", no_argument, NULL, '5' },
    { "L2_ICM", no_argument, NULL, '6' },
    { "L2_DCM", no_argument, NULL, '7' },
    { "TLB_IM", no_argument, NULL, '8' },
    { "TLB_DM", no_argument, NULL, '9' },
#endif
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
        printf("-N N          Use N blocks per matrix (default N = %llu)\n", N);
        printf("--align N     Align memory buffer on N byte boundary (default N = %llu)\n", alignment);
        printf("--loops N     Repeat each access test N times (default N = %llu)\n", loops);
#ifdef HAVE_PAPI
        printf("--TOT_INS     Measure total instructions\n");
        printf("--TOT_CYC     Measure total cycles\n");
        printf("--RES_STL     Measure stalled cycles\n");
        printf("--L1_ICM      Measure L1 instruction misses\n");
        printf("--L1_DCM      Measure L1 data misses\n");
        printf("--L2_ICM      Measure L2 instruction misses\n");
        printf("--L2_DCM      Measure L2 data misses\n");
        printf("--TLB_IM      Measure TLB instruction misses\n");
        printf("--TLB_DM      Measure TLB data misses\n");
#endif
        return 0;
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'a':
        alignment = strtol(optarg, NULL, 10);
        break;

      case 'l':
        loops = strtol(optarg, NULL, 10);
        break;

#ifdef HAVE_PAPI
      case '1':
        load_TOT_INS = 1;
        break;

      case '2':
        load_TOT_CYC = 1;
        break;

      case '3':
        load_RES_STL = 1;
        break;

      case '4':
        load_L1_ICM = 1;
        break;

      case '5':
        load_L1_DCM = 1;
        break;

      case '6':
        load_L2_ICM = 1;
        break;

      case '7':
        load_L2_DCM = 1;
        break;

      case '8':
        load_TLB_IM = 1;
        break;

      case '9':
        load_TLB_DM = 1;
        break;
#endif

      default:
        printf("unknown command line argument\n");
        return -1;
        break;
    }
  }

#ifdef HAVE_PAPI
  /* Do some PAPI. */
  if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
  {
    printf("can not initialize PAPI\n");
    exit(1);
  }

  papi_values = (long long*) malloc(sizeof(long long)*9);
  if ((papi_errorcode = PAPI_set_granularity(PAPI_GRN_MIN) != PAPI_OK)) { print_papi_error("failed to set granularity",  papi_errorcode); }
  if ((papi_errorcode = PAPI_create_eventset(&papi_events) != PAPI_OK)) { print_papi_error("failed to create eventset",  papi_errorcode); }
  if ((papi_errorcode = PAPI_assign_eventset_component(papi_events, 0) != PAPI_OK)) { print_papi_error("failed to assign component to eventset",  papi_errorcode); }

  if (load_TOT_INS && (papi_errorcode = PAPI_add_event(papi_events, PAPI_TOT_INS))) { print_papi_error("failed to add PAPI_TOT_INS", papi_errorcode); }
  if (load_TOT_CYC && (papi_errorcode = PAPI_add_event(papi_events, PAPI_TOT_CYC))) { print_papi_error("failed to add PAPI_TOT_CYC", papi_errorcode); }
  if (load_RES_STL && (papi_errorcode = PAPI_add_event(papi_events, PAPI_RES_STL))) { print_papi_error("failed to add PAPI_RES_STL", papi_errorcode); }

  if (load_L1_ICM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L1_ICM))) { print_papi_error("failed to add PAPI_L1_ICM",  papi_errorcode); }
  if (load_L1_DCM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L1_DCM))) { print_papi_error("failed to add PAPI_L1_DCM",  papi_errorcode); }
  if (load_L2_ICM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L2_ICM))) { print_papi_error("failed to add PAPI_L2_ICM",  papi_errorcode); }
  if (load_L2_DCM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L2_DCM))) { print_papi_error("failed to add PAPI_L2_DCM",  papi_errorcode); }
  if (load_TLB_IM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_TLB_IM))) { print_papi_error("failed to add PAPI_TLB_IM",  papi_errorcode); }
  if (load_TLB_DM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_TLB_DM))) { print_papi_error("failed to add PAPI_TLB_DM",  papi_errorcode); }
#endif

#ifdef HAVE_POSIX_MEMALIGN
  if ((result = posix_memalign((void**) &A, alignment, sizeof(float)*16*N)) != 0)
  {
    switch (result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("There was insufficient memory to fulfill the allocation request.\n");
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }

  if ((result = posix_memalign((void**) &B, alignment, sizeof(float)*16*N)) != 0)
  {
    switch (result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("There was insufficient memory to fulfill the allocation request.\n");
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }

  if ((result = posix_memalign((void**) &C, alignment, sizeof(float)*16*N)) != 0)
  {
    switch (result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("There was insufficient memory to fulfill the allocation request.\n");
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }

  if ((result = posix_memalign((void**) &stream, alignment, sizeof(float)*3*16*N)) != 0)
  {
    switch (result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("There was insufficient memory to fulfill the allocation request.\n");
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }
#else
  if ((A = (float*) malloc(sizeof(float)*16*N)) == NULL)
  {
    printf("error allocating A\n");
    exit(1);
  }

  if ((B = (float*) malloc(sizeof(float)*16*N)) == NULL)
  {
    printf("error allocating B\n");
    exit(1);
  }

  if ((C = (float*) malloc(sizeof(float)*16*N)) == NULL)
  {
    printf("error allocating C\n");
    exit(1);
  }

  if ((stream = (float*) malloc(sizeof(float)*3*16*N)) == NULL)
  {
    printf("error allocating stream\n");
    exit(1);
  }
#endif

  /* Load random data into stream. */
  for (index = 0; index < N; index++)
  {
    A_block = &A[index*16];
    B_block = &B[index*16];
    C_block = &C[index*16];

    //A_block = &stream[index*3*16];
    //B_block = &stream[index*3*16+16];
    //C_block = &stream[index*3*16+32];

    for (i = 0; i < N_BLOCK; i++) {
      for (j = 0; j < N_BLOCK; j++)
      {
        A_block[i*N_BLOCK+j] = rand()/(float) RAND_MAX;
        B_block[i*N_BLOCK+j] = rand()/(float) RAND_MAX;
        C_block[i*N_BLOCK+j] = rand()/(float) RAND_MAX;
      }
    }
  }

  printf("stream has %llu matrix blocks\n", N);
  if (3*16*N*sizeof(float) < 1024) printf("stream is %llu bytes\n", 3*16*N*sizeof(float));
  else if (3*16*N*sizeof(float) < 1024*1024) printf("stream is %1.1f kB\n", 3*16*N*sizeof(float)/1024.);
  else if (3*16*N*sizeof(float) < 1024*1024*1024) printf("stream is %1.1f MB\n", 3*16*N*sizeof(float)/1024./1024.);
  else printf("stream is %1.1f GB\n", 3*16*N*sizeof(float)/1024./1024./1024.);
  printf("stream uses %llu pages of 4 kB\n", 3*16*N*sizeof(float)/4096);
  printf("if this were a matrix multiply, the underlying matrices would have dimension %f x %f\n", pow(N, 1/3.)*N_BLOCK, pow(N, 1/3.)*N_BLOCK);
  printf("loops = %llu\n", loops);

  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);

#ifdef HAVE_PAPI
  papi_errorcode = PAPI_num_events(papi_events);
  if (papi_errorcode == 0)
  {
    printf("[PAPI] no counters specified\n");
  }

  else if (papi_errorcode < 0)
  {
    print_papi_error("failed to count number of counters in eventset", papi_errorcode);
  }

  else
  {
    if ((papi_errorcode = PAPI_start(papi_events))) { print_papi_error("failed start", papi_errorcode); }
  }
#endif

  for (loop = 0; loop < loops; loop++)
  {
    for (index = 0; index < N; index++)
    {
      /* Load pointer to blocks. */
      A_block = &A[index*16];
      B_block = &B[index*16];
      C_block = &C[index*16];

      /* Do some prefetching. */
      _mm_prefetch((void*) &A[(index+6)*16], _MM_HINT_T0);
      _mm_prefetch((void*) &B[(index+6)*16], _MM_HINT_T0);
      _mm_prefetch((void*) &C[(index+6)*16], _MM_HINT_T0);

      //A_block = &stream[index*3*16];
      //B_block = &stream[index*3*16+16];
      //C_block = &stream[index*3*16+32];

      if (A_block[0] == 0) { break; }
      if (B_block[0] == 0) { break; }
      if (C_block[0] == 0) { break; }
    }
  }

#ifdef HAVE_PAPI
  if (PAPI_num_events(papi_events) > 0)
  {
    if (PAPI_stop(papi_events, papi_values) != PAPI_OK) { printf("failed to stop\n"); exit(1); }
    papi_value_index = 0;

    if (load_TOT_INS)
    {
      printf("[PAPI] total instructions = %lli, per block = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) N);
      papi_value_index++;
    }

    if (load_TOT_CYC)
    {
      printf("[PAPI] total cycles = %lli, per block = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) N);
      papi_value_index++;
    }

    if (load_RES_STL)
    {
      printf("[PAPI] cycles stalled = %lli, per block = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) N);
      papi_value_index++;
    }

    if (load_L1_ICM)
    {
      printf("[PAPI] instruction L1 misses = %lli, per block = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) N);
      papi_value_index++;
    }

    if (load_L1_DCM)
    {
      printf("[PAPI] data L1 misses = %lli, per block = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) N);
      papi_value_index++;
    }

    if (load_L2_ICM)
    {
      printf("[PAPI] instruction L2 misses = %lli, per block = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) N);
      papi_value_index++;
    }

    if (load_L2_DCM)
    {
      printf("[PAPI] data L2 misses = %lli, per block = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) N);
      papi_value_index++;
    }

    if (load_TLB_IM)
    {
      printf("[PAPI] total instruction TLB misses = %lli, per block = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) N);
      papi_value_index++;
    }

    if (load_TLB_DM)
    {
      printf("[PAPI] total data TLB misses = %lli, per block = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) N);
      papi_value_index++;
    }
  }

  PAPI_cleanup_eventset(papi_events);
  PAPI_destroy_eventset(&papi_events);
  free(papi_values);
#endif

  getrusage(RUSAGE_SELF, &rusage_stop);
  gettimeofday(&stop, NULL);
  walltime = (stop.tv_sec-start.tv_sec+(stop.tv_usec-start.tv_usec)/1.0e6)/loops;
  usertime = ((rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6)/loops;
  systime = ((rusage_stop.ru_stime.tv_sec-rusage_start.ru_stime.tv_sec)+(rusage_stop.ru_stime.tv_usec-rusage_start.ru_stime.tv_usec)/(double) 1e6)/loops;
  flops = (pow(N, 1/3.)*N_BLOCK)*(pow(N, 1/3.)*N_BLOCK)*(2.*pow(N, 1/3.)*N_BLOCK+1.)/usertime;
  printf("performance: total walltime %f s, usertime %f s, systime %f s, ", walltime*loops, usertime*loops, systime*loops);
  printf("per loop walltime %e s, usertime %e s, systime %e, ", walltime, usertime, systime);
  if (flops < 1000*1000*1000)
  {
    printf("%1.2f Mflop/s\n", flops/1000./1000.);
  }

  else
  {
    printf("%1.2f Gflop/s\n", flops/1000./1000./1000.);
  }

  free(A);
  free(B);
  free(C);
}
