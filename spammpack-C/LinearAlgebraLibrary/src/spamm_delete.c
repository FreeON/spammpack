#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Delete a matrix.
 *
 * @param A The matrix to delete.
 */
void
spamm_delete (struct spamm_t *A)
{
  assert(A != NULL);

  /* Recurse down and free() nodes. */
  LOG2_INFO("deleting matrix\n");
  if (A->root != NULL)
  {
    spamm_delete_node(&A->root);
  }

  A->M = 0;
  A->N = 0;

  A->M_padded = 0;
  A->N_padded = 0;

  A->M_child = 0;
  A->N_child = 0;

  A->M_block = 0;
  A->N_block = 0;
}
