#include <spamm.h>
#include <stdlib.h>

int
main ()
{
  int i, j;
  struct spamm_t A;

  spamm_new(10, 10, 1, 1, 2, 2, 0.0, &A);

  for (i = 0; i < 10; ++i) {
    for (j = 0; j < 10; ++j)
    {
      //spamm_set(i, j, rand()/(float) RAND_MAX, &A);
      spamm_set(i, j, spamm_dense_index(i, j, 10, 10), &A);
    }
  }

  //spamm_print_tree(&A);
  spamm_print_spamm(&A);

  spamm_tree_pack(0, 100, i_mask, &A);

  return 0;
}
