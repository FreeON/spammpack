#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

void
spamm_set_element (const int i, const int j, const double Aij, struct spamm_node_t *node)
{
  assert(node != NULL);

  int l, k;

  //spamm_log("node [%i-%i[ x [%i-%i[\n", __FILE__, __LINE__, node->M_lower, node->M_upper, node->N_lower, node->N_upper);
  //spamm_print_node(node);

  /* Check if we are at the block level. */
  if ((node->M_upper-node->M_lower) == node->M_block && (node->N_upper-node->N_lower) == node->N_block)
  {
    //spamm_log("reached block level\n", __FILE__, __LINE__);
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

    //spamm_log("storing A[%i][%i] in block_dense[%i][%i]\n", __FILE__, __LINE__, i, j, i-node->M_lower, j-node->N_lower);
    node->block_dense[spamm_dense_index(i-node->M_lower, j-node->N_lower, node->N_block)] = Aij;
  }

  else
  {
    if (node->child[0] == NULL)
    {
      /* Create children nodes. */
      //spamm_log("creating children\n", __FILE__, __LINE__);
      for (l = 0; l < 2; ++l) {
        for (k = 0; k < 2; ++k)
        {
          spamm_new_node(&(node->child[spamm_dense_index(l, k, 2)]));
          node->child[spamm_dense_index(l, k, 2)]->M_lower = node->M_lower+l*(node->M_upper-node->M_lower)/2;
          node->child[spamm_dense_index(l, k, 2)]->M_upper = node->M_lower+(l+1)*(node->M_upper-node->M_lower)/2;
          node->child[spamm_dense_index(l, k, 2)]->N_lower = node->N_lower+k*(node->N_upper-node->N_lower)/2;
          node->child[spamm_dense_index(l, k, 2)]->N_upper = node->N_lower+(k+1)*(node->N_upper-node->N_lower)/2;

          node->child[spamm_dense_index(l, k, 2)]->M_block = node->M_block;
          node->child[spamm_dense_index(l, k, 2)]->N_block = node->N_block;
        }
      }
    }

    /* Recurse. */
    if (i < (node->M_lower+(node->M_upper-node->M_lower)/2) && j < (node->N_lower+(node->N_upper-node->N_lower)/2))
    {
      return spamm_set_element(i, j, Aij, node->child[spamm_dense_index(0, 0, 2)]);
    }

    else if (i < (node->M_lower+(node->M_upper-node->M_lower)/2))
    {
      return spamm_set_element(i, j, Aij, node->child[spamm_dense_index(0, 1, 2)]);
    }

    else if (j < (node->N_lower+(node->N_upper-node->N_lower)/2))
    {
      return spamm_set_element(i, j, Aij, node->child[spamm_dense_index(1, 0, 2)]);
    }

    else
    {
      return spamm_set_element(i, j, Aij, node->child[spamm_dense_index(1, 1, 2)]);
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
    //spamm_log("creating new root node\n", __FILE__, __LINE__);
    spamm_new_node(&(A->root));
    A->root->M_lower = 0;
    A->root->M_upper = A->M_padded;
    A->root->N_lower = 0;
    A->root->N_upper = A->N_padded;
    A->root->M_block = A->M_block;
    A->root->N_block = A->N_block;
  }

  //spamm_log("recursing to root node to store A[%i][%i] = %f\n", __FILE__, __LINE__, i, j, Aij);
  spamm_set_element(i, j, Aij, A->root);
}
