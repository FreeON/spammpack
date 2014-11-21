#include <spamm.h>

int
main ()
{
  int result = 0;

  const unsigned int N = 1000;
  const unsigned int chunk_tier = 5;
  const short use_linear_tree = 1;

  struct spamm_matrix_t *A;

  A = spamm_new(N, N, chunk_tier, use_linear_tree, row_major);

  spamm_expand(A);

  spamm_delete(&A);

  return result;
}
