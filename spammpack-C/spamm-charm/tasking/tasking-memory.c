#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define COLUMN_MAJOR(i, j, N) ((i)+(N)*(j))

void work (const int N,
           const int P,
           const int tier,
           const int depth) {

}

int main (int argc, char **argv) {

  int N = 256;
  int P = 1;
  int depth = 0;

  double *A;

  const char *short_options = "h";
  char c;

  while((c = getopt(argc, argv, short_options)) != -1) {
    switch(c) {
    case 'h':
      printf("Usage:\n");
      printf("\n");
      printf("-h        This help\n");
      printf("-N N      The matrix size (currently %d)\n", N);
      printf("-P N      The fan-out per tier (currently %d)\n", P);
      printf("-d depth  The recursion depth (currently %d)\n", depth);
      exit(0);
      break;

    case 'N':
      N = strtoll(optarg, NULL, 10);
      break;

    case 'd':
      depth = strtol(optarg, NULL, 10);
      break;

    case 'P':
      P = strtol(optarg, NULL, 10);
      break;

    default:
      fprintf(stderr, "illegal argument\n");
      exit(1);
      break;
    }
  }

  A = calloc(N*N, sizeof(double));
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      A[COLUMN_MAJOR(i, j, N)] = rand()/(double) RAND_MAX;
    }
  }

#pragma omp parallel
  {
#pragma omp single
    {
#pragma omp task untied shared(N)
      {
        work(N, P, 0, depth);
      }
#pragma omp taskwait
    }
  }

  return 0;
}
