#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

struct linear_t
{
  unsigned int index;
  floating_point_t *block_dense;
};

int
main ()
{
  int i, j;
  floating_point_t A_dense[16] = {
    1, 5, 9,  13,
    2, 6, 10, 14,
    3, 7, 11, 15,
    4, 8, 12, 16
  };

  struct spamm_t A;
  struct linear_t *linear_A;

  printf("A =\n");
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < 4; ++j)
    {
      printf(" %2.0f", A_dense[i+j*4]);
    }
    printf("\n");
  }

  /* Construct quadtree. */
  spamm_new(4, 4, 1, 1, 2, 2, 0.0, &A);
  spamm_dense_to_spamm(4, 4, 1, 1, 2, 2, 0.0, A_dense, &A);
  linear_A = (struct linear_t*) malloc(sizeof(struct linear_t)*A.number_nonzero_blocks);

  /* Print tree. */
  spamm_print_tree(&A);

  /* Free quadtree. */
  free(linear_A);
  spamm_delete(&A);

  return 0;
}
