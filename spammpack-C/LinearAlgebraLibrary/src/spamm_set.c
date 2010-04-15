#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

void
spamm_set_element (const int i, const int j, const double Aij, struct spamm_node_t *node)
{
  assert(node != NULL);

  int l, k;

  /* Check if we are at the block level. */
  if ((node->M_upper-node->M_lower) == node->M_block && (node->N_upper-node->N_lower) == node->N_block)
  {
    if (node->block_dense == NULL)
    {
      node->block_dense = (double*) malloc(sizeof(double)*node->M_block*node->N_block);
      for (l = 0; l < node->M_block; ++l) {
        for (k = 0; k < node->N_block; ++k)
        {
          node->block_dense[spamm_dense_index(l, k, node->N_block)] = 0.0;
        }
      }
    }

    node->block_dense[spamm_dense_index(i-node->M_lower, j-node->N_lower, node->N_block)] = Aij;
  }

  else
  {
    if (node->child == NULL)
    {
      /* Allocate children nodes. */
      node->child = (struct spamm_node_t**) malloc(sizeof(struct spamm_node_t*)*node->M_child*node->N_child);

      /* Create children nodes. */
      for (l = 0; l < node->M_child; ++l) {
        for (k = 0; k < node->N_child; ++k)
        {
          spamm_new_node(&(node->child[spamm_dense_index(l, k, node->N_child)]));
          node->child[spamm_dense_index(l, k, node->N_child)]->M_lower = node->M_lower+l*(node->M_upper-node->M_lower)/node->M_child;
          node->child[spamm_dense_index(l, k, node->N_child)]->M_upper = node->M_lower+(l+1)*(node->M_upper-node->M_lower)/node->M_child;
          node->child[spamm_dense_index(l, k, node->N_child)]->N_lower = node->N_lower+k*(node->N_upper-node->N_lower)/node->N_child;
          node->child[spamm_dense_index(l, k, node->N_child)]->N_upper = node->N_lower+(k+1)*(node->N_upper-node->N_lower)/node->N_child;

          node->child[spamm_dense_index(l, k, node->N_child)]->M_child = node->M_child;
          node->child[spamm_dense_index(l, k, node->N_child)]->N_child = node->N_child;

          node->child[spamm_dense_index(l, k, node->N_child)]->M_block = node->M_block;
          node->child[spamm_dense_index(l, k, node->N_child)]->N_block = node->N_block;
        }
      }
    }

    /* Recurse. */
    for (l = 0; l < node->M_child; ++l) {
      for (k = 0; k < node->N_child; ++k)
      {
        if (i >= (node->M_lower+(node->M_upper-node->M_lower)*l/node->M_child) &&
            i < (node->M_lower+(node->M_upper-node->M_lower)*(l+1)/node->M_child) &&
            j >= (node->N_lower+(node->N_upper-node->N_lower)*k/node->N_child) &&
            j < (node->N_lower+(node->N_upper-node->N_lower)*(k+1)/node->N_child))
        {
          return spamm_set_element(i, j, Aij, node->child[spamm_dense_index(l, k, node->N_child)]);
        }
      }
    }
  }
}

void
spamm_set (const int i, const int j, const double Aij, struct spamm_t *A)
{
  assert(A != NULL);

  if (i < 0)     { spamm_log("i < 0\n",  __FILE__, __LINE__); exit(1); }
  if (j < 0)     { spamm_log("j < 0\n",  __FILE__, __LINE__); exit(1); }
  if (i >= A->M) { spamm_log("i >= M\n", __FILE__, __LINE__); exit(1); }
  if (i >= A->N) { spamm_log("i >= N\n", __FILE__, __LINE__); exit(1); }

  /* Recursively find the leaf node that stores this element. */
  if (A->root == NULL)
  {
    spamm_new_node(&(A->root));
    A->root->M_lower = 0;
    A->root->M_upper = A->M_padded;
    A->root->N_lower = 0;
    A->root->N_upper = A->N_padded;

    A->root->M_child = A->M_child;
    A->root->N_child = A->N_child;

    A->root->M_block = A->M_block;
    A->root->N_block = A->N_block;
  }

  spamm_set_element(i, j, Aij, A->root);
}
