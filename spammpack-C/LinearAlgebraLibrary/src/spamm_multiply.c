#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/* Computes
 *
 * C_node = alpha*A_node*B_node + beta*C_node
 */
void
spamm_multiply_node (const double alpha, const struct spamm_node_t *A_node, const struct spamm_node_t *B_node, const double beta, struct spamm_node_t **C_node)
{
  int i, j, k;

  for (i = 0; i < (*C_node)->M_child; ++i) {
    for (j = 0; j < (*C_node)->N_child; ++j) {
      for (k = 0; k < A_node->N_child; ++k)
      {
      }
    }
  }
}

/* Computes
 *
 * C = alpha*A*B + beta*C
 */
void
spamm_multiply (const double alpha, const struct spamm_t *A, const struct spamm_t *B, const double beta, struct spamm_t *C)
{
  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

  if (A->N != B->M)
  {
    spamm_log("matrix size mismatch, A->N = %i, B->M = %i\n", __FILE__, __LINE__, A->N, B->M);
    exit(1);
  }

  if (A->M != C->M)
  {
    spamm_log("matrix size mismatch, A->M = %i, C->M = %i\n", __FILE__, __LINE__, A->M, C->M);
    exit(1);
  }

  if (B->N != C->N)
  {
    spamm_log("matrix size mismatch, B->N = %i, C->N = %i\n", __FILE__, __LINE__, B->N, C->N);
    exit(1);
  }

  if (A->N_child != B->M_child)
  {
    spamm_log("matrix child size mismatch, A->N_child = %i, B_child->M = %i\n", __FILE__, __LINE__, A->N_child, B->M_child);
    exit(1);
  }

  if (A->M_child != C->M_child)
  {
    spamm_log("matrix child size mismatch, A->M_child = %i, C->M_child = %i\n", __FILE__, __LINE__, A->M_child, C->M_child);
    exit(1);
  }

  if (B->N_child != C->N_child)
  {
    spamm_log("matrix child size mismatch, B->N_child = %i, C->N_child = %i\n", __FILE__, __LINE__, B->N_child, C->N_child);
    exit(1);
  }

  if (A->N_block != B->M_block)
  {
    spamm_log("matrix block size mismatch, A->N_block = %i, B_block->M = %i\n", __FILE__, __LINE__, A->N_block, B->M_block);
    exit(1);
  }

  if (A->M_block != C->M_block)
  {
    spamm_log("matrix block size mismatch, A->M_block = %i, C->M_block = %i\n", __FILE__, __LINE__, A->M_block, C->M_block);
    exit(1);
  }

  if (B->N_block != C->N_block)
  {
    spamm_log("matrix block size mismatch, B->N_block = %i, C->N_block = %i\n", __FILE__, __LINE__, B->N_block, C->N_block);
    exit(1);
  }

  if (alpha == 0.0 || A->root == NULL || B->root == NULL)
  {
    /* Only multiply C by beta. */
    spamm_log("[FIXME]", __FILE__, __LINE__);
    exit(1);
  }

  else
  {
    spamm_multiply_node(alpha, A->root, B->root, beta, &(C->root));
  }
}
