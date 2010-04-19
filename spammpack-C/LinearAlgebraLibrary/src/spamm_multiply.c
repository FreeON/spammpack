#include "spamm.h"
#include "config.h"
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

  if (*C_node == NULL)
  {
    /* Create new node. */
    spamm_new_node(C_node);
    (*C_node)->M_lower = A_node->M_lower;
    (*C_node)->M_upper = A_node->M_upper;
    (*C_node)->N_lower = B_node->N_lower;
    (*C_node)->N_upper = B_node->N_upper;

    (*C_node)->M_child = A_node->M_child;
    (*C_node)->N_child = B_node->N_child;

    (*C_node)->threshold = A_node->threshold;

    (*C_node)->M_block = A_node->M_block;
    (*C_node)->N_block = B_node->N_block;
  }

  if (A_node->child != NULL && B_node->child != NULL)
  {
    if ((*C_node)->child == NULL)
    {
      (*C_node)->child = (struct spamm_node_t**) malloc(sizeof(struct spamm_node_t*)*(*C_node)->M_child*(*C_node)->N_child);
      for (i = 0; i < (*C_node)->M_child; ++i) {
        for (j = 0; j < (*C_node)->N_child; ++j)
        {
          spamm_new_node(&((*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]));
          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->M_lower = (*C_node)->M_lower+i*((*C_node)->M_upper-(*C_node)->M_lower)/(*C_node)->M_child;
          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->M_upper = (*C_node)->M_lower+(i+1)*((*C_node)->M_upper-(*C_node)->M_lower)/(*C_node)->M_child;
          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->N_lower = (*C_node)->N_lower+j*((*C_node)->N_upper-(*C_node)->N_lower)/(*C_node)->N_child;
          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->N_upper = (*C_node)->N_lower+(j+1)*((*C_node)->N_upper-(*C_node)->N_lower)/(*C_node)->N_child;

          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->M_child = (*C_node)->M_child;
          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->N_child = (*C_node)->N_child;

          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->threshold = (*C_node)->threshold;

          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->M_block = (*C_node)->M_block;
          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->N_block = (*C_node)->N_block;
        }
      }
    }

    for (i = 0; i < (*C_node)->M_child; ++i) {
      for (j = 0; j < (*C_node)->N_child; ++j) {
        for (k = 0; k < A_node->N_child; ++k)
        {
          spamm_multiply_node(alpha, A_node->child[spamm_dense_index(i, k, A_node->M_child, A_node->N_child)],
              B_node->child[spamm_dense_index(k, j, B_node->M_child, B_node->N_child)], beta,
              &(*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]);
        }
      }
    }
  }

  if (A_node->block_dense != NULL && B_node->block_dense != NULL)
  {
    if ((*C_node)->block_dense == NULL)
    {
      /* Create empty dense block. */
      (*C_node)->block_dense = (double*) malloc(sizeof(double)*(*C_node)->M_block*(*C_node)->N_block);
      for (i = 0; i < (*C_node)->M_block; ++i) {
        for (j = 0; j < (*C_node)->N_block; ++j)
        {
          (*C_node)->block_dense[spamm_dense_index(i, j, (*C_node)->M_block, (*C_node)->N_block)] = 0;
        }
      }
    }

#ifdef DGEMM
    DGEMM("N", "N", &(A_node->M_block), &(B_node->N_block), &(A_node->N_block),
        &alpha, A_node->block_dense, &(A_node->M_block), B_node->block_dense, &(B_node->M_block),
        &beta, (*C_node)->block_dense, &((*C_node)->M_block));
#else
    for (i = 0; i < (*C_node)->M_block; ++i) {
      for (j = 0; j < (*C_node)->N_block; ++j) {
        for (k = 0; k < A_node->M_block; ++k)
        {
          (*C_node)->block_dense[spamm_dense_index(i, j, (*C_node)->M_block, (*C_node)->N_block)]
            = alpha*A_node->block_dense[spamm_dense_index(i, k, A_node->M_block, A_node->N_block)]*B_node->block_dense[spamm_dense_index(k, j, B_node->M_block, B_node->N_block)]
            + beta*(*C_node)->block_dense[spamm_dense_index(i, j, (*C_node)->M_block, (*C_node)->N_block)];
        }
      }
    }
#endif
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

  if (C->root != NULL)
  {
    spamm_log("[FIXME] can not handle pre-existing C\n", __FILE__, __LINE__);
    exit(1);
  }

  if (beta != 1.0)
  {
    spamm_log("[FIXME] can not handle (beta = %e) != 1.0\n", __FILE__, __LINE__, beta);
    exit(1);
  }

  if ((A->root == NULL || B->root == NULL) && C->root == NULL)
  {
    /* Nothing to be done. */
    return;
  }

  else if ((A->root == NULL || B->root == NULL) && C->root != NULL)
  {
    spamm_add_node(0.0, NULL, beta, &(C->root));
  }

  spamm_multiply_node(alpha, A->root, B->root, beta, &(C->root));
}
