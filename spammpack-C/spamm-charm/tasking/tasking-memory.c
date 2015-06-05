#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#ifdef INTEL_TASK_API
#include <ittnotify.h>
__itt_domain *domain;
#endif

#define COLUMN_MAJOR(i, j, N) ((i)+(N)*(j))

/** Recursively multiply two matrices.
 *
 * @param N The matrix size.
 */
void work (const int N, const int i_lower, const int i_upper,
           const int j_lower, const int j_upper, const int k_lower,
           const int k_upper,
           const double *A,
           const double *B,
           double *C,
           const int tier,
           const int depth)
{
  if(tier == depth) {
#ifdef INTEL_TASK_API
    char task_name[100];
    sprintf(task_name, "leaf[%d:%d,%d:%d,%d:%d]", i_lower, i_upper, j_lower, j_upper,
            k_lower, k_upper);
    __itt_string_handle *leaf_task_handle = __itt_string_handle_create(task_name);
    __itt_task_begin(domain, __itt_null, __itt_null, leaf_task_handle);
#endif
    /* Block multiply. */
    int width = i_upper-i_lower;
    int A_offset = width*width*COLUMN_MAJOR(i_lower/width, k_lower/width, N/width);
    int B_offset = width*width*COLUMN_MAJOR(k_lower/width, j_lower/width, N/width);
    int C_offset = width*width*COLUMN_MAJOR(i_lower/width, j_lower/width, N/width);
    /* printf("[%d:%d,%d:%d,%d:%d]\n", i_lower, i_upper, j_lower, j_upper, k_lower, k_upper); */
    /* printf("A_offset = %d\n", A_offset); */
    /* printf("B_offset = %d\n", B_offset); */
    /* printf("C_offset = %d\n", C_offset); */
    for(int i = 0; i < width; i++) {
      for(int j = 0; j < width; j++) {
        for(int k = 0; k < width; k++) {
          C[C_offset+COLUMN_MAJOR(i, j, width)] +=
            A[A_offset+COLUMN_MAJOR(i, k, width)]*B[B_offset+COLUMN_MAJOR(k, j, width)];
        }
      }
    }
#ifdef INTEL_TASK_API
    __itt_task_end(domain);
#endif
  } else {
    /* Descend. */
    int half_width = (i_upper-i_lower)/2;
    for(int k = 0; k < 2; k++) {
      for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
#ifdef INTEL_TASK_API
            char task_name[100];
            sprintf(task_name, "tree[%d:%d,%d:%d,%d:%d]", i_lower, i_upper, j_lower, j_upper,
                    k_lower, k_upper);
            __itt_string_handle *tree_handle = __itt_string_handle_create(task_name);
#endif
#pragma omp task untied firstprivate(i, j, k)
          {
#ifdef INTEL_TASK_API
            __itt_task_begin(domain, __itt_null, __itt_null, tree_handle);
#endif
            work(N,
                 i_lower+i*half_width, i_lower+(i+1)*half_width,
                 j_lower+j*half_width, j_lower+(j+1)*half_width,
                 k_lower+k*half_width, k_lower+(k+1)*half_width,
                 A, B, C, tier+1, depth);
#ifdef INTEL_TASK_API
            __itt_task_end(domain);
#endif
          }
        }
      }
#pragma omp taskwait
    }
  }
}

int main (int argc, char **argv)
{
  int N = 256;
  int N_basic = 256;
  int depth = 0;
  int repeat = 1;

  short depth_set = 0;
  short N_basic_set = 0;

  struct timespec timer_resolution;
  struct timespec start_time;
  struct timespec end_time;

  const char *short_options = "hN:d:b:r:";
  char c;

  while((c = getopt(argc, argv, short_options)) != -1) {
    switch(c) {
    case 'h':
      printf("Usage:\n");
      printf("\n");
      printf("-h        This help\n");
      printf("-N N      The matrix size (currently %d)\n", N);
      printf("-d depth  The recursion depth (currently %d)\n", depth);
      printf("-b N      The basic matrix size (currently %d)\n", N_basic);
      printf("-r R      Number of repetitions (currently %d)\n", repeat);
      exit(0);
      break;

    case 'N':
      N = strtoll(optarg, NULL, 10);
      break;

    case 'd':
      depth = strtol(optarg, NULL, 10);
      depth_set = 1;
      break;

    case 'b':
      N_basic = strtol(optarg, NULL, 10);
      N_basic_set = 1;
      break;

    case 'r':
      repeat = strtol(optarg, NULL, 10);
      break;

    default:
      exit(1);
      break;
    }
  }

  if(depth_set) {
    if(N_basic_set) {
      printf("Set only depth or N_basic, not both\n");
      exit(1);
    }
    N_basic = N;
    for(int tier = 0; tier < depth; tier++) {
      N_basic >>= 1;
    }
    printf("N_basic = %d\n", N_basic);
  }

  if(N_basic_set) {
    if(depth_set) {
      printf("Set only depth or N_basic, not both\n");
      exit(1);
    }
    int N_temp = N;
    for(depth = 0; N_temp > N_basic; depth++) {
      N_temp >>= 1;
    }
    if(N_temp != N_basic) {
      printf("logic error, N_basic incorrect, found %d\n", N_temp);
      exit(1);
    }
  }

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
  for(int tier = 1; tier <= depth; tier++) {
    number_tree_tasks += number_leaf_tasks;
    number_leaf_tasks *= 8;
  }
  printf("creating %d tree and %d leaf tasks, respectively\n",
         number_tree_tasks, number_leaf_tasks);
  printf("depth = %d, N_basic = %d\n", depth, N_basic);
  printf("initializing random %dx%d matrix... ", N, N);
  double *A = calloc(N*N, sizeof(double));
  double *C = calloc(N*N, sizeof(double));
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      A[COLUMN_MAJOR(i, j, N)] = rand()/(double) RAND_MAX;
    }
  }
  printf("done\n");

#ifdef INTEL_TASK_API
  domain = __itt_domain_create("tasking domain");
#endif

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
#ifdef _OPENMP
      printf("running with %d threads\n", omp_get_num_threads());
#else
      printf("running serial version\n");
#endif
      for(int i = 0; i < repeat; i++) {
#pragma omp task untied
        {
#ifdef INTEL_TASK_API
          __itt_string_handle *top_tree_handle = __itt_string_handle_create("top tree");
          __itt_task_begin(domain, __itt_null, __itt_null, top_tree_handle);
#endif
          work(N, 0, N, 0, N, 0, N, A, A, C, 0, depth);
#ifdef INTEL_TASK_API
          __itt_task_end(domain);
#endif
        }
#pragma omp taskwait
      }
    }
  }

  if(clock_gettime(CLOCK_MONOTONIC, &end_time) != 0) {
    fprintf(stderr, "can not get end time\n");
    exit(1);
  }

  printf("timer resolution: %e ns\n",
         (double) timer_resolution.tv_nsec
         + (double) timer_resolution.tv_sec*1e9);
  printf("elapsed time: %e s\n",
         (double) (end_time.tv_nsec-start_time.tv_nsec) * 1e-9
         + (double) (end_time.tv_sec-start_time.tv_sec));

  return 0;
}
