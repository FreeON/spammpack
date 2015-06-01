#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>
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
           const int depth) {

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
    for(int i = 0; i < width; i++) {
      for(int j = 0; j < width; j++) {
        for(int k = 0; k < width; k++) {
          C[COLUMN_MAJOR(i, j, width)] +=
            A[COLUMN_MAJOR(i, k, width)]*B[COLUMN_MAJOR(k, j, width)];
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
#pragma omp task firstprivate(i, j, k)
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

int main (int argc, char **argv) {

  int N = 256;
  int depth = 0;

  const char *short_options = "hN:d:";
  char c;

  while((c = getopt(argc, argv, short_options)) != -1) {
    switch(c) {
    case 'h':
      printf("Usage:\n");
      printf("\n");
      printf("-h        This help\n");
      printf("-N N      The matrix size (currently %d)\n", N);
      printf("-d depth  The recursion depth (currently %d)\n", depth);
      exit(0);
      break;

    case 'N':
      N = strtoll(optarg, NULL, 10);
      break;

    case 'd':
      depth = strtol(optarg, NULL, 10);
      break;

    default:
      exit(1);
      break;
    }
  }

  int number_leaf_tasks = 1;
  int number_tree_tasks = 0;
  int N_basic = N;
  for(int tier = 1; tier <= depth; tier++) {
    number_tree_tasks += number_leaf_tasks;
    number_leaf_tasks *= 8;
    N_basic /= 2;
  }
  printf("creating %d tree and %d leaf tasks, respectively\n",
         number_tree_tasks, number_leaf_tasks);
  printf("N_basic = %d\n", N_basic);
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

#pragma omp parallel
  {
#pragma omp single
    {
#ifdef _OPENMP
      printf("running with %d threads\n", omp_get_num_threads());
#else
      printf("running serial version\n");
#endif
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

  return 0;
}
