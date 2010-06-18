#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** \private Calculate
 *
 * \f$B_{\mathrm node} = \alpha A_{\mathrm node} + \beta B_{\mathrm node}\f$
 *
 * This is the recursive part. Use spamm_add() instead.
 *
 * @param alpha The factor \f$\alpha\f$.
 * @param A_node Matrix node \f$A\f$.
 * @param beta The factor \f$\beta\f$.
 * @param B_node Matrix node \f$B\f$.
 */
void
spamm_add_node (const float_t alpha, const struct spamm_node_t *A_node, const float_t beta, struct spamm_node_t **B_node)
{
  int i, j;
  struct spamm_node_t *new_node;
  struct spamm_ll_iterator_t *iterator_A, *iterator_B;
  struct spamm_ll_node_t *linear_node_A, *linear_node_B;
  struct spamm_linear_quadtree_t *linear_A, *linear_B;

  /* We distinguish between the following cases:
   *
   * (1) neither A nor B exist --> We are done.
   * (2) A exists but B doesn't --> We need to create a new node in B.
   */

  if (A_node == NULL && *B_node == NULL)
  {
    /* We are done here. */
    return;
  }

  if (A_node != NULL && *B_node == NULL)
  {
    /* We need to add to B. */
    *B_node = spamm_new_node();

    (*B_node)->tier = A_node->tier;
    (*B_node)->tree_depth = A_node->tree_depth;

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

  /* Decide on how to recurse further. */
  if (A_node != NULL && A_node->child != NULL && (*B_node)->child == NULL)
  {
    /* We need to recurse further. A still has children nodes, but B doesn't.
     * We therefore create children nodes in B.
     */
    (*B_node)->child = (struct spamm_node_t**) malloc(sizeof(struct spamm_node_t*)*(*B_node)->M_child*(*B_node)->N_child);
    for (i = 0; i < A_node->M_child; ++i) {
      for (j = 0; j < A_node->N_child; ++j)
      {
        (*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)] = spamm_new_node();
        new_node = (*B_node)->child[spamm_dense_index(i, j, (*B_node)->M_child, (*B_node)->N_child)];

        new_node->tier = (*B_node)->tier+1;
        new_node->tree_depth = (*B_node)->tree_depth;

        new_node->M_lower = (*B_node)->M_lower+i*((*B_node)->M_upper-(*B_node)->M_lower)/(*B_node)->M_child;
        new_node->M_upper = (*B_node)->M_lower+(i+1)*((*B_node)->M_upper-(*B_node)->M_lower)/(*B_node)->M_child;
        new_node->N_lower = (*B_node)->N_lower+j*((*B_node)->N_upper-(*B_node)->N_lower)/(*B_node)->N_child;
        new_node->N_upper = (*B_node)->N_lower+(j+1)*((*B_node)->N_upper-(*B_node)->N_lower)/(*B_node)->N_child;

        new_node->M_child = (*B_node)->M_child;
        new_node->N_child = (*B_node)->N_child;

        new_node->threshold = (*B_node)->threshold;

        new_node->linear_tier = (*B_node)->linear_tier;

        new_node->M_block = (*B_node)->M_block;
        new_node->N_block = (*B_node)->N_block;
      }
    }
  }

  if (A_node != NULL && A_node->block_dense != NULL && (*B_node)->block_dense == NULL)
  {
    /* We can stop recursing. A has a dense block but B doesn't. We therefore
     * create an empty dense block in B. */
    (*B_node)->block_dense = (float_t*) malloc(sizeof(float_t)*(*B_node)->M_block*(*B_node)->N_block);
    for (i = 0; i < (*B_node)->M_block; ++i) {
      for (j = 0; j < (*B_node)->N_block; ++j)
      {
        (*B_node)->block_dense[spamm_dense_index(i, j, (*B_node)->M_block, (*B_node)->N_block)] = 0;
      }
    }
  }

  if (A_node != NULL && A_node->linear_quadtree != NULL && (*B_node)->linear_quadtree == NULL)
  {
    /* Create a new linear quadtree in B. */
    (*B_node)->linear_quadtree = spamm_ll_new();
    (*B_node)->linear_quadtree_memory = spamm_mm_new(A_node->linear_quadtree_memory->chunksize);
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

    /* A doesn't exist, only multiply B with beta. */
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

    /* A doesn't exist. Only multiply B with beta. */
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

  else if ((*B_node)->linear_quadtree != NULL)
  {
    /* Multiply B with beta. */
    iterator_B = spamm_ll_iterator_new((*B_node)->linear_quadtree);
    for (linear_node_B = spamm_ll_iterator_first(iterator_B); linear_node_B != NULL; linear_node_B = spamm_ll_iterator_next(iterator_B))
    {
      linear_B = linear_node_B->data;
      for (i = 0; i < (*B_node)->M_block; ++i) {
        for (j = 0; j < (*B_node)->N_block; ++j)
        {
          linear_B->block_dense[spamm_dense_index(i, j, (*B_node)->M_block, (*B_node)->N_block)] *= beta;
        }
      }
    }
    spamm_ll_iterator_delete(&iterator_B);

    /* Add to B to A. */
    if (A_node != NULL && A_node->linear_quadtree != NULL)
    {
      iterator_A = spamm_ll_iterator_new(A_node->linear_quadtree);
      for (linear_node_A = spamm_ll_iterator_first(iterator_A); linear_node_A != NULL; linear_node_A = spamm_ll_iterator_next(iterator_A))
      {
        linear_A = linear_node_A->data;

        /* Search B to find whether we have a matching block that we can add.
         * If not, we need to create a new one.
         */
        linear_B = NULL;
        iterator_B = spamm_ll_iterator_new((*B_node)->linear_quadtree);
        for (linear_node_B = spamm_ll_iterator_first(iterator_B); linear_node_B != NULL; linear_node_B = spamm_ll_iterator_next(iterator_B))
        {
          linear_B = linear_node_B->data;
          if (linear_A->index == linear_B->index)
          {
            /* Found matching block in B. */
            break;
          }
        }
        spamm_ll_iterator_delete(&iterator_B);

        if (linear_B == NULL || (linear_B != NULL && linear_A->index != linear_B->index))
        {
          /* We didn't find a matching block in B. Create a new one. */
          linear_B = (struct spamm_linear_quadtree_t*) spamm_mm_allocate(sizeof(struct spamm_linear_quadtree_t)+(*B_node)->M_block*(*B_node)->N_block*sizeof(float_t)+8, (*B_node)->linear_quadtree_memory);
          spamm_ll_append(linear_B, (*B_node)->linear_quadtree);
          linear_B->block_dense = (float_t*) (((void*) linear_B)+sizeof(struct spamm_linear_quadtree_t));
          linear_B->index = linear_A->index;
          for (i = 0; i < (*B_node)->M_block; ++i) {
            for (j = 0; j < (*B_node)->N_block; ++j)
            {
              linear_B->block_dense[spamm_dense_index(i, j, (*B_node)->M_block, (*B_node)->N_block)] = 0.0;
            }
          }
        }

        /* Add blocks of A and B. */
        for (i = 0; i < (*B_node)->M_block; ++i) {
          for (j = 0; j < (*B_node)->N_block; ++j)
          {
            linear_B->block_dense[spamm_dense_index(i, j, (*B_node)->M_block, (*B_node)->N_block)] += alpha*linear_A->block_dense[spamm_dense_index(i, j, A_node->M_block, A_node->N_block)];
          }
        }
      }
      spamm_ll_iterator_delete(&iterator_A);
    }
  }
}

/** Calculate
 *
 * \f$B = \alpha A + \beta B\f$
 *
 * @param alpha The factor \f$\alpha\f$.
 * @param A Matrix \f$A\f$.
 * @param beta The factor \f$\beta\f$.
 * @param B Matrix \f$B\f$.
 */
void
spamm_add (const float_t alpha, const struct spamm_t *A, const float_t beta, struct spamm_t *B)
{
  assert(A != NULL);
  assert(B != NULL);

  if (A->M != B->M)
  {
    LOG_FATAL("matrix size mismatch, A->M = %i, B->M = %i\n", A->M, B->M);
    exit(1);
  }

  if (A->N != B->N)
  {
    LOG_FATAL("matrix size mismatch, A->N = %i, B->N = %i\n", A->N, B->N);
    exit(1);
  }

  if (A->M_child != B->M_child)
  {
    LOG_FATAL("matrix child size mismatch, A->M_child = %i, B->M_child = %i\n", A->M_child, B->M_child);
    exit(1);
  }

  if (A->N_child != B->N_child)
  {
    LOG_FATAL("matrix child size mismatch, A->N_child = %i, B->N_child = %i\n", A->N_child, B->N_child);
    exit(1);
  }

  if (A->M_block != B->M_block)
  {
    LOG_FATAL("matrix block size mismatch, A->M_block = %i, B->M_block = %i\n", A->M_block, B->M_block);
    exit(1);
  }

  if (A->N_block != B->N_block)
  {
    LOG_FATAL("matrix block size mismatch, A->N_block = %i, B->N_block = %i\n", A->N_block, B->N_block);
    exit(1);
  }

  if (A->linear_tier != B->linear_tier)
  {
    LOG_FATAL("linear_tier mismatch, A->linear_tier = %u, B->linear_tier = %u\n", A->linear_tier, B->linear_tier);
    exit(1);
  }

  spamm_add_node(alpha, A->root, beta, &(B->root));
}
