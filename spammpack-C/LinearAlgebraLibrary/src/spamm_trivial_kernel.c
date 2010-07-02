#include "spamm.h"

void
spamm_sgemm_trivial (const floating_point_t alpha,
    const struct spamm_node_t *A_node,
    const struct spamm_node_t *B_node,
    struct spamm_node_t *C_node)
{
  unsigned int i, j, k;

  for (i = 0; i < C_node->M_block; i++) {
    for (j = 0; j < C_node->N_block; j++) {
      for (k = 0; k < A_node->N_block; k++)
      {
        C_node->block_dense[spamm_dense_index(i, j, C_node->M_block, C_node->N_block)]
          += alpha*A_node->block_dense[spamm_dense_index(i, k, A_node->M_block, A_node->N_block)]
          *B_node->block_dense[spamm_dense_index(k, j, B_node->M_block, B_node->N_block)];
        //C_node->block_dense[i+j*C_node->M_block]
        //  += alpha*A_node->block_dense[i+k*A_node->M_block]
        //  *B_node->block_dense[k+j*B_node->M_block];
      }
    }
  }
}
