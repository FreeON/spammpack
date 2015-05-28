#ifdef _OPENMP
#include <omp.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

int work (const long int N, const int P, const int tier, const int depth) {

  int result;
  long int i, j;

  result = 0;

  if(tier < depth) {
    for(i = 0; i < P; i++) {
#pragma omp task untied default(shared) firstprivate(result)
      {
        result = work(N, P, tier+1, depth);
      }
    }
#pragma omp taskwait
  }
  else {
    result = 0;
    for(i = 0; i < N; i++) {
      j = i;
      while(j > 10) j -= 10;
      result += j;
      while(result > 15) result -= 15;
    }
  }

  return result;
}

int main (int argc, char **argv) {

  long int N = 100000;
  int P = 1;
  int depth = 0;

  long int number_tasks = 1;

  struct timespec timer_resolution;
  struct timespec start_time;
  struct timespec end_time;

  int i;
  int result;
  char c;
  char *short_options = "hN:d:P:";

  while((c = getopt(argc, argv, short_options)) != -1) {
    switch(c) {
    case 'N':
      N = strtol(optarg, NULL, 10);
      break;

    case 'd':
      depth = strtol(optarg, NULL, 10);
      break;

    case 'P':
      P = strtol(optarg, NULL, 10);
      break;

    case 'h':
      printf("Usage:\n");
      printf("\n");
      printf("-h        This help\n");
      printf("-N N      The size of the tasks\n");
      printf("-P N      The fan-out per tier\n");
      printf("-d depth  The recursion depth\n");
      exit(0);
      break;

    default:
      fprintf(stderr, "illegal argument\n");
      exit(-1);
    }
  }

  for(i = 1; i <= depth; i++) {
    number_tasks *= P;
  }

  printf("fan-out: %d, depth: %d, %ld work tasks\n", P, depth, number_tasks);

  if(clock_getres(CLOCK_MONOTONIC, &timer_resolution) != 0) {
    fprintf(stderr, "can not get timer resolution\n");
    exit(1);
  }
  if(clock_gettime(CLOCK_MONOTONIC, &start_time) != 0) {
    fprintf(stderr, "can not get start time\n");
    exit(1);
  }
#pragma omp parallel
  {
#pragma omp single
    {
#pragma omp task untied default(shared)
      {
        result = work(N, P, 0, depth);
      }
#pragma omp taskwait
    }
  }
  if(clock_gettime(CLOCK_MONOTONIC, &end_time) != 0) {
    fprintf(stderr, "can not get end time\n");
    exit(1);
  }

  if(result < 0) {
    printf("result is < 0\n");
  }

  printf("timer resolution: %e ns\n", (double) timer_resolution.tv_nsec
         + (double) timer_resolution.tv_sec*1e9);
  printf("elapsed time: %e s\n", (double) (end_time.tv_nsec-start_time.tv_nsec) * 1e-9
         + (double) (end_time.tv_sec-start_time.tv_sec));

  return 0;
}
