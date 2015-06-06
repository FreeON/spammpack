#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "tasking-memory.h"

#define COLUMN_MAJOR(i, j, N) ((i)+(N)*(j))

/** Recursively multiply two matrices.
 *
 * @param N The matrix size.
 */
void work(const int N, const int i_lower, const int i_upper,
          const int j_lower, const int j_upper, const int k_lower,
          const int k_upper,
          const double *A,
          const double *B,
          double *C,
          const int tier,
          const int depth)
{
  if (tier == depth) {
    intel_task_begin("leaf[%d:%d,%d:%d,%d:%d]", i_lower, i_upper, j_lower, j_upper,
                     k_lower, k_upper);
    /* Block multiply. */
    int width = i_upper-i_lower;
    int A_offset = width*width*COLUMN_MAJOR(i_lower/width, k_lower/width, N/width);
    int B_offset = width*width*COLUMN_MAJOR(k_lower/width, j_lower/width, N/width);
    int C_offset = width*width*COLUMN_MAJOR(i_lower/width, j_lower/width, N/width);
    for (int i = 0; i < width; i++) {
      for (int j = 0; j < width; j++) {
        for (int k = 0; k < width; k++) {
          C[C_offset+COLUMN_MAJOR(i, j, width)] +=
            A[A_offset+COLUMN_MAJOR(i, k, width)]*B[B_offset+COLUMN_MAJOR(k, j, width)];
        }
      }
    }
    intel_task_end();
  } else {
    /* Descend. */
    int half_width = (i_upper-i_lower)/2;
    for (int k = 0; k < 2; k++) {
      for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
#pragma omp task untied firstprivate(i, j, k)
          {
            intel_task_begin("tree[%d:%d,%d:%d,%d:%d]", i_lower, i_upper, j_lower, j_upper,
                             k_lower, k_upper);
            work(N,
                 i_lower+i*half_width, i_lower+(i+1)*half_width,
                 j_lower+j*half_width, j_lower+(j+1)*half_width,
                 k_lower+k*half_width, k_lower+(k+1)*half_width,
                 A, B, C, tier+1, depth);
            intel_task_end();
          }
        }
      }
#pragma omp taskwait
    }
  }
}

int main(int argc, char **argv)
{
  int N = 256;
  int N_chunk = 256;
  int N_basic = 256;
  int depth = 0;
  int repeat = 1;
  int tile_repeat = 1;

  short depth_set = 0;
  short N_chunk_set = 0;
  short N_basic_set = 0;
  short repeat_set = 0;

  struct timespec timer_resolution;
  struct timespec start_time;
  struct timespec end_time;

  const char *short_options = "hN:c:b:d:r:";
  char c;

  while ((c = getopt(argc, argv, short_options)) != -1) {
    switch (c) {
    case 'h':
      printf("Usage:\n");
      printf("\n");
      printf("-h        This help\n");
      printf("-N N      The matrix size (currently %d)\n", N);
      printf("-c N      The chunk matrix size (currently %d)\n", N_chunk);
      printf("-b N      The basic matrix size (currently %d)\n", N_basic);
      printf("-d depth  The recursion depth (currently %d)\n", depth);
      printf("-r R      Number of repetitions (currently %d)\n", repeat);
      exit(0);
      break;

    case 'N':
      N = strtoll(optarg, NULL, 10);
      break;

    case 'c':
      N_chunk = strtol(optarg, NULL, 10);
      N_chunk_set = 1;
      break;

    case 'b':
      N_basic = strtol(optarg, NULL, 10);
      N_basic_set = 1;
      break;

    case 'd':
      depth = strtol(optarg, NULL, 10);
      depth_set = 1;
      break;

    case 'r':
      repeat = strtol(optarg, NULL, 10);
      repeat_set = 1;
      break;

    default:
      exit(1);
      break;
    }
  }

  if (N_chunk_set) {
    int N_temp = N_chunk;
    while (N_temp < N) {
      N_temp <<= 1;
    }
    if (N_temp > N) {
      printf("logic error, N_chunk = %d incorrect\n", N_chunk);
      exit(1);
    }
  } else {
    N_chunk = N;
  }

  if (!(depth_set || N_basic_set)) {
    N_basic = N_chunk;
  }

  if (depth_set) {
    if (N_basic_set) {
      printf("Set only depth or N_basic, not both\n");
      exit(1);
    }
    N_basic = N_chunk;
    for (int tier = 0; tier < depth; tier++) {
      N_basic >>= 1;
    }
    printf("N_basic = %d\n", N_basic);
  }

  if (N_basic_set) {
    if (depth_set) {
      printf("Set only depth or N_basic, not both\n");
      exit(1);
    }
    int N_temp = N_chunk;
    for (depth = 0; N_temp > N_basic; depth++) {
      N_temp >>= 1;
    }
    if (N_temp != N_basic) {
      printf("logic error, N_basic = %d incorrect\n", N_basic);
      exit(1);
    }
  }

  if (N_basic > N_chunk) {
    printf("N_basic > N_chunk\n");
    exit(1);
  }
  if (N_chunk > N) {
    printf("N_chunk > N\n");
    exit(1);
  }

  printf("N = %d, N_chunk = %d, N_basic = %d, depth = %d\n", N, N_chunk, N_basic, depth);

  if (N%N_chunk != 0) {
    printf("logic error, N%%N_chunk = %d\n", N%N_chunk);
    exit(1);
  }
  tile_repeat = (N/N_chunk)*(N/N_chunk)*(N/N_chunk);
  printf("setting tile_repeat to %d\n", tile_repeat);

  /* Startup threads. */
#pragma omp parallel
  {
#pragma omp single
    {
      printf("starting...\n");
    }
  }

  int number_leaf_tasks = 1;
  int number_tree_tasks = 0;
  for (int tier = 1; tier <= depth; tier++) {
    number_tree_tasks += number_leaf_tasks;
    number_leaf_tasks *= 8;
  }
  printf("creating %d tree and %d leaf tasks, respectively\n",
         number_tree_tasks, number_leaf_tasks);
  printf("initializing random %dx%d matrix... ", N_chunk, N_chunk);
  double *A = calloc(N_chunk*N_chunk, sizeof(double));
  double *C = calloc(N_chunk*N_chunk, sizeof(double));
  for (int i = 0; i < N_chunk; i++) {
    for (int j = 0; j < N_chunk; j++) {
      A[COLUMN_MAJOR(i, j, N_chunk)] = rand()/(double) RAND_MAX;
    }
  }
  printf("done\n");

  intel_task_domain_create();

  if (clock_getres(CLOCK_MONOTONIC, &timer_resolution) != 0) {
    fprintf(stderr, "can not get timer resolution\n");
    exit(1);
  }
  if (clock_gettime(CLOCK_MONOTONIC, &start_time) != 0) {
    fprintf(stderr, "can not get start time\n");
    exit(1);
  }

#pragma omp parallel
  {
#pragma omp single
    {
#ifdef _OPENMP
      printf("running with %d threads\n", omp_get_num_threads());
#else
      printf("running serial version\n");
#endif
      for (int i = 0; i < repeat; i++ ) {
        for (int tile = 0; tile < tile_repeat; tile++) {
#pragma omp task untied
          {
            intel_task_begin("top tree");
            work(N_chunk, 0, N_chunk, 0, N_chunk, 0, N_chunk, A, A, C, 0, depth);
            intel_task_end();
          }
#pragma omp taskwait
        }
      }
    }
  }

  if (clock_gettime(CLOCK_MONOTONIC, &end_time) != 0) {
    fprintf(stderr, "can not get end time\n");
    exit(1);
  }

  printf("timer resolution: %e ns\n",
        (double) timer_resolution.tv_nsec
        + (double) timer_resolution.tv_sec*1e9);
  printf("elapsed time: %e s\n",
        ((double) (end_time.tv_nsec-start_time.tv_nsec) * 1e-9
        + (double) (end_time.tv_sec-start_time.tv_sec))/repeat);

  return 0;
}
