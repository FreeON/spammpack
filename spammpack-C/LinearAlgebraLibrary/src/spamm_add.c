/** @file */

#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** \private Calculate
 *
 * B_node = alpha*A_node + beta*B_node
 *
 * This is the recursive part. Use spamm_add() instead.
 *
 * @param alpha The factor alpha.
 * @param A_node Matrix node A.
 * @param beta The factor beta.
 * @param B_node Matrix node B.
 */
void
spamm_add_node (const float_t alpha, const struct spamm_node_t *A_node, const float_t beta, struct spamm_node_t **B_node)
{
  int i, j;

  if (A_node == NULL && *B_node == NULL)
  {
    /* We are done here. */
    return;
  }

  if (A_node != NULL && *B_node == NULL)
  {
    /* We need to add to B. */
    spamm_new_node(B_node);

    (*B_node)->tier = A_node->tier;

    (*B_node)->M_lower = A_node->M_lower;
    (*B_node)->M_upper = A_node->M_upper;
    (*B_node)->N_lower = A_node->N_lower;
    (*B_node)->N_upper = A_node->N_upper;

    (*B_node)->M_child = A_node->M_child;
    (*B_node)->N_child = A_node->N_child;

    (*B_node)->threshold = A_node->threshold;

    (*B_node)->linear_tier = A_node->linear_tier;

    (*B_node)->M_block = A_node->M_block;
    (*B_node)->N_block = A_node->N_block;
  }

  if (A_node != NULL && A_node->child != NULL && (*B_node)->child == NULL)
  {
    (*B_node)->child = (struct spamm_node_t**) malloc(sizeof(struct spamm_node_t*)*(*B_node)->M_child*(*B_node)->N_child);
    for (i = 0; i < A_node->M_child; ++i) {
      for (j = 0; j < A_node->N_child; ++j)
      {
        spamm_new_node(&((*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]));

        (*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]->tier = (*B_node)->tier+1;

        (*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]->M_lower = (*B_node)->M_lower+i*((*B_node)->M_upper-(*B_node)->M_lower)/(*B_node)->M_child;
        (*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]->M_upper = (*B_node)->M_lower+(i+1)*((*B_node)->M_upper-(*B_node)->M_lower)/(*B_node)->M_child;
        (*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]->N_lower = (*B_node)->N_lower+j*((*B_node)->N_upper-(*B_node)->N_lower)/(*B_node)->N_child;
        (*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]->N_upper = (*B_node)->N_lower+(j+1)*((*B_node)->N_upper-(*B_node)->N_lower)/(*B_node)->N_child;

        (*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]->M_child = (*B_node)->M_child;
        (*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]->N_child = (*B_node)->N_child;

        (*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]->threshold = (*B_node)->threshold;

        (*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]->linear_tier = (*B_node)->linear_tier;

        (*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]->M_block = (*B_node)->M_block;
        (*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]->N_block = (*B_node)->N_block;
      }
    }
  }

  if (A_node != NULL && A_node->block_dense != NULL && (*B_node)->block_dense == NULL)
  {
    /* Create empty dense block. */
    (*B_node)->block_dense = (float_t*) malloc(sizeof(float_t)*(*B_node)->M_block*(*B_node)->N_block);
    for (i = 0; i < (*B_node)->M_block; ++i) {
      for (j = 0; j < (*B_node)->N_block; ++j)
      {
        (*B_node)->block_dense[spamm_dense_index(i, j, (*B_node)->M_block, (*B_node)->N_block)] = 0;
      }
    }
  }

  if ((*B_node)->child != NULL)
  {
    /* Recurse further down in A & B. */
    if (A_node != NULL && A_node->child != NULL)
    {
      for (i = 0; i < A_node->M_child; ++i) {
        for (j = 0; j < A_node->N_child; ++j)
        {
          spamm_add_node(alpha, A_node->child[spamm_dense_index(i, j, A_node->M_child, A_node->N_child)], beta, &((*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]));
        }
      }
    }

    else
    {
      for (i = 0; i < (*B_node)->M_child; ++i) {
        for (j = 0; j < (*B_node)->N_child; ++j)
        {
          spamm_add_node(alpha, NULL, beta, &((*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)]));
        }
      }
    }
  }

  else if ((*B_node)->block_dense != NULL)
  {
    /* Add dense blocks. */
    if (A_node != NULL && A_node->block_dense != NULL)
    {
      for (i = 0; i < (*B_node)->M_block; ++i) {
        for (j = 0; j < (*B_node)->N_block; ++j)
        {
          (*B_node)->block_dense[spamm_dense_index(i, j, (*B_node)->M_block, (*B_node)->N_block)] = alpha*A_node->block_dense[spamm_dense_index(i, j, A_node->M_block, A_node->N_block)]
            + beta*(*B_node)->block_dense[spamm_dense_index(i, j, (*B_node)->M_block, (*B_node)->N_block)];
        }
      }
    }

    else
    {
      for (i = 0; i < (*B_node)->M_block; ++i) {
        for (j = 0; j < (*B_node)->N_block; ++j)
        {
          (*B_node)->block_dense[spamm_dense_index(i, j, (*B_node)->M_block, (*B_node)->N_block)] = beta*(*B_node)->block_dense[spamm_dense_index(i, j, (*B_node)->M_block, (*B_node)->N_block)];
        }
      }
    }
  }

  else
  {
    //spamm_log("else?\n", __FILE__, __LINE__);
  }
}

/** Calculate
 *
 * B = alpha*A + beta*B
 *
 *
 * @param alpha The factor alpha.
 * @param A Matrix A.
 * @param beta The factor beta.
 * @param B Matrix B.
 */
void
spamm_add (const float_t alpha, const struct spamm_t *A, const float_t beta, struct spamm_t *B)
{
  assert(A != NULL);
  assert(B != NULL);

  if (A->M != B->M)
  {
    spamm_log("matrix size mismatch, A->M = %i, B->M = %i\n", __FILE__, __LINE__, A->M, B->M);
    exit(1);
  }

  if (A->N != B->N)
  {
    spamm_log("matrix size mismatch, A->N = %i, B->N = %i\n", __FILE__, __LINE__, A->N, B->N);
    exit(1);
  }

  if (A->M_child != B->M_child)
  {
    spamm_log("matrix child size mismatch, A->M_child = %i, B->M_child = %i\n", __FILE__, __LINE__, A->M_child, B->M_child);
    exit(1);
  }

  if (A->N_child != B->N_child)
  {
    spamm_log("matrix child size mismatch, A->N_child = %i, B->N_child = %i\n", __FILE__, __LINE__, A->N_child, B->N_child);
    exit(1);
  }

  if (A->M_block != B->M_block)
  {
    spamm_log("matrix block size mismatch, A->M_block = %i, B->M_block = %i\n", __FILE__, __LINE__, A->M_block, B->M_block);
    exit(1);
  }

  if (A->N_block != B->N_block)
  {
    spamm_log("matrix block size mismatch, A->N_block = %i, B->N_block = %i\n", __FILE__, __LINE__, A->N_block, B->N_block);
    exit(1);
  }

  spamm_add_node(alpha, A->root, beta, &(B->root));
}
