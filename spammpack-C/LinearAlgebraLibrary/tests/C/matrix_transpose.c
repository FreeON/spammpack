#include <lal.h>

int
main ()
{
  int M = 10;
  int N = 20;

  int i, j;

  lal_matrix_t *A;
  lal_matrix_t *A_transpose;

  lal_allocate(M, N, &A);

  lal_rand(A);

  A_transpose = lal_transpose(A);

  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      if (lal_get(i, j, A) != lal_get(j, i, A_transpose))
      {
        return -1;
      }
    }
  }

  return 0;
}
