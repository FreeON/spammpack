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
  floating_point_t A_dense[64] = {
    1, 9,  17, 25, 33, 41, 49, 57,
    2, 10, 18, 26, 34, 42, 50, 58,
    3, 11, 19, 27, 35, 43, 51, 59,
    4, 12, 20, 28, 36, 44, 52, 60,
    5, 13, 21, 29, 37, 45, 53, 61,
    6, 14, 22, 30, 38, 46, 54, 62,
    7, 15, 23, 31, 39, 47, 55, 63,
    8, 16, 24, 32, 40, 48, 56, 64
  };

  struct spamm_t A;
  struct linear_t *linear_A;

  printf("A =\n");
  for (i = 0; i < 8; ++i) {
    for (j = 0; j < 8; ++j)
    {
      printf(" %2.0f", A_dense[i+j*8]);
    }
    printf("\n");
  }

  /* Construct quadtree. */
  spamm_new(8, 8, &A);
  spamm_dense_to_spamm(8, 8, 'N', A_dense, &A);
  linear_A = (struct linear_t*) malloc(sizeof(struct linear_t)*A.number_nonzero_blocks);

  /* Print tree. */
  spamm_print_tree(&A);

  /* Print SpAMM. */
  spamm_print_spamm(&A);

  /* Free quadtree. */
  free(linear_A);
  spamm_delete(&A);

  return 0;
}
