#include "spamm.h"

/** Create the tree necessary to store an element at A(i,j).
 *
 * @param i The row index.
 * @param j The column index.
 * @param A The matrix.
 */
void
spamm_create_tree (const unsigned int i, const unsigned int j, struct spamm_t *A)
{
  assert(A != NULL);

  if (i >= A->M) { LOG_FATAL("(i = %i) >= (M = %i)\n", i, A->M); exit(1); }
  if (j >= A->N) { LOG_FATAL("(j = %i) >= (N = %i)\n", j, A->N); exit(1); }

  /* Recursively find the leaf node that stores this element. */
  if (A->root == NULL)
  {
    A->root = spamm_new_childnode(0, A->tree_depth,
        0, A->N_padded, 0, A->N_padded,
        0, 0, 0, 0, A->kernel_tier, NULL, NULL);
  }
}
