#include <spamm.h>
#include <math.h>
#include <stdlib.h>

int
main ()
{
  int result = 0;

  unsigned int N = 1200;
  unsigned int i, j;
  floating_point_t *A_dense;
  floating_point_t norm, norm2;
  struct spamm_t A;

  A_dense = (floating_point_t*) malloc(sizeof(floating_point_t)*N*N);
  norm2 = 0;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      A_dense[spamm_dense_index(i, j, N, N)] = rand()/(floating_point_t) RAND_MAX;
      norm2 += A_dense[spamm_dense_index(i, j, N, N)]*A_dense[spamm_dense_index(i, j, N, N)];
    }
  }
  norm = sqrt(norm2);

  spamm_new(N, N, &A);
  spamm_dense_to_spamm(N, N, A_dense, &A);

  if (norm != A.norm)
  {
    printf("norm = %f, A.norm = %f, diff = %e, rel. diff = %e\n",
        norm, A.norm, fabs(norm-A.norm), fabs(norm-A.norm)/norm);
    result = -1;
  }

  return result;
}
