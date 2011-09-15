#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "spamm.h"
#include "spamm_naive.h"

/** Create a new naive matrix object.
 *
 * @param M Number of rows of dense input matrix.
 * @param N Number of columns of dense input matrix.
 * @param blocksize The size of the dense matrix blocks.
 *
 * @return A pointer to the matrix.
 */
struct spamm_naive_t *
spamm_naive_new (const unsigned int M, const unsigned int N, const unsigned int blocksize)
{
  struct spamm_naive_t *A = NULL;
  double x, x_M, x_N;

  if (M <= 0)
  {
    fprintf(stderr, "M <= 0\n");
    exit(1);
  }

  if (N <= 0)
  {
    fprintf(stderr, "N <= 0\n");
    exit(1);
  }

  /* Allocate memory. */
  A = calloc(1, sizeof(struct spamm_naive_t));

  /* Store the blocksize. */
  A->blocksize = blocksize;

  /* Pad to powers of M_child x N_child. */
  x_M = (log(M) > log(blocksize) ? log(M) - log(blocksize) : 0)/log(2);
  x_N = (log(N) > log(blocksize) ? log(N) - log(blocksize) : 0)/log(2);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  /* The ceil() function can lead to a depth that is one tier too large
   * because of numerical errors in the calculation of x. We need to check
   * whether the depth is appropriate.
   */
  A->depth = (unsigned int) ceil(x);

  /* Double check depth. */
  if (A->depth >= 1 && ((int) (blocksize*pow(2, A->depth-1)) >= M && (int) (blocksize*pow(2, A->depth-1)) >= N))
  {
    (A->depth)--;
  }

  /* Set matrix size. */
  A->M = M;
  A->N = N;

  /* Set padded matrix size. */
  A->N_padded = (int) (blocksize*pow(2, A->depth));

  return A;
}

/** Allocate a new node of a naive matrix tree.
 *
 * @param tier The tier this node will be on.
 *
 * @return A pointer to the newly allocated node.
 */
struct spamm_naive_node_t *
spamm_naive_new_node (const unsigned int tier)
{
  struct spamm_naive_node_t *node = NULL;

  node = calloc(1, sizeof(struct spamm_naive_node_t));

  return node;
}

/** @brief Get an element from a naive matrix.
 *
 * If the matri is NULL, this function returns 0.
 *
 * @param i The row index.
 * @param j The column index.
 * @param A The matrix.
 *
 * @return The matrix element Aij.
 */
double
spamm_naive_get (const unsigned int i, const unsigned int j, const struct spamm_naive_t *A)
{
  struct spamm_naive_node_t **node = NULL;

  if (A == NULL)
  {
    return 0;
  }

  node = (struct spamm_naive_node_t**) &(A->root);

  while (1)
  {
    if (*node == NULL)
    {
      return 0;
    }

    if ((*node)->M_upper-(*node)->M_lower == A->blocksize)
    {
      /* Get the matrix element. */
      if ((*node)->data == NULL) { return 0.0; }
      else
      {
        return (*node)->data[spamm_index_column_major(i-(*node)->M_lower, j-(*node)->N_lower, A->blocksize, A->blocksize)];
      }
    }

    else
    {
      if (i < (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
          j < (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
      {
        node = &((*node)->child[0]);
      }

      else if (i <  (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
               j >= (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
      {
        node = &((*node)->child[1]);
      }

      else if (i >= (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
               j <  (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
      {
        node = &((*node)->child[2]);
      }

      else if (i >= (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
               j >= (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
      {
        node = &((*node)->child[3]);
      }
    }
  }
}

/** Recursively set a matrix element.
 *
 * @param i The row index.
 * @param j The column index.
 * @param Aij The value of the matrix element A(i,j).
 * @param tier The tier the node is on.
 * @param node The node.
 */
void
spamm_naive_set_recursive (const unsigned int i, const unsigned int j, const float Aij,
    const unsigned int M_lower,
    const unsigned int M_upper,
    const unsigned int N_lower,
    const unsigned int N_upper,
    const int blocksize,
    const int tier,
    struct spamm_naive_node_t **node)
{
  if (*node == NULL)
  {
    /* Allocate new node. */
    *node = spamm_naive_new_node(tier+1);

    (*node)->M_lower = M_lower;
    (*node)->M_upper = M_upper;
    (*node)->N_lower = N_lower;
    (*node)->N_upper = N_upper;

    (*node)->tier = tier;
  }

  if ((*node)->M_upper-(*node)->M_lower == blocksize)
  {
    /* Store the matrix element. */
    if ((*node)->data == NULL)
    {
      (*node)->blocksize = blocksize;
      (*node)->data = calloc(blocksize*blocksize, sizeof(float));
    }

    /* sgemm() loves column major. */
    (*node)->data[spamm_index_column_major(i-(*node)->M_lower, j-(*node)->N_lower, blocksize, blocksize)] = Aij;

    /* Update norm. */
    (*node)->norm2 += Aij*Aij;
    (*node)->norm   = sqrt((*node)->norm2);
  }

  else
  {
    if (i < (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
        j < (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
    {
      spamm_naive_set_recursive(i, j, Aij,
          (*node)->M_lower, (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2,
          (*node)->N_lower, (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2,
          blocksize, tier+1, &((*node)->child[0]));
    }

    else if (i <  (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
        j >= (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
    {
      spamm_naive_set_recursive(i, j, Aij,
          (*node)->M_lower, (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2,
          (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2, (*node)->N_upper,
          blocksize, tier+1, &((*node)->child[1]));
    }

    else if (i >= (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
        j <  (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
    {
      spamm_naive_set_recursive(i, j, Aij,
          (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2, (*node)->M_upper,
          (*node)->N_lower, (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2,
          blocksize, tier+1, &((*node)->child[2]));
    }

    else if (i >= (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2 &&
        j >= (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2)
    {
      spamm_naive_set_recursive(i, j, Aij,
          (*node)->M_lower+((*node)->M_upper-(*node)->M_lower)/2, (*node)->M_upper,
          (*node)->N_lower+((*node)->N_upper-(*node)->N_lower)/2, (*node)->N_upper,
          blocksize, tier+1, &((*node)->child[3]));
    }

    /* Update norm. */
    (*node)->norm2 += Aij*Aij;
    (*node)->norm   = sqrt((*node)->norm2);
  }
}

/** Set an element in a naive matrix.
 *
 * @param i The row index.
 * @param j The column index.
 * @param Aij The value of the matrix element A(i,j).
 * @param A The matrix.
 */
void
spamm_naive_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_naive_t *A)
{
  if (A == NULL)
  {
    printf("[%s:%i] A can not be NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  /* Recurse the tree to find the dense data block to store the matrix element
   * in.
   */
  if (Aij == 0.0) { return; }

  /* Recursively set the matrix element. */
  spamm_naive_set_recursive(i, j, Aij, 0, A->N_padded, 0, A->N_padded, A->blocksize, 0, &(A->root));
}

/** Multiply to matrices, i.e. \f$ C = \alpha A \times B + \beta C\f$.
 *
 * @param tolerance The SpAMM tolerance of this product.
 * @param alpha The paramater \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param beta The paramater \f$\beta\f$.
 * @param C The matrix \f$C\f$.
 * @param timer The timer to use.
 * @param sgemm The external sgemm function to use.
 */
void
spamm_naive_multiply (const float tolerance,
    const float alpha, struct spamm_naive_t *A, struct spamm_naive_t *B,
    const float beta, struct spamm_naive_t *C,
    struct spamm_timer_t *timer,
    void (*sgemm) (char *, char *, int *, int *, int *, float *, float *, int *, float *, int *, float *, float *, int *))
{
  if (A == NULL)
  {
    printf("[%s:%i] A is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (B == NULL)
  {
    printf("[%s:%i] B is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (C == NULL)
  {
    printf("[%s:%i] C is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (A->blocksize != B->blocksize)
  {
    printf("[%s:%i] A->blocksize != B->blocksize\n", __FILE__, __LINE__);
    exit(1);
  }

  if (A->blocksize != C->blocksize)
  {
    printf("[%s:%i] A->blocksize != C->blocksize\n", __FILE__, __LINE__);
    exit(1);
  }

  /* Multiply C by beta. */
  spamm_naive_multiply_scalar(beta, C->root);

  /* Multiply A and B. */
  if (A->root != NULL && B->root != NULL && C->root == NULL)
  {
    C->root = spamm_naive_new_node(0);
    C->root->M_lower = 0;
    C->root->M_upper = A->root->M_upper;
    C->root->N_lower = 0;
    C->root->N_upper = A->root->N_upper;
  }

  spamm_naive_multiply_matrix(tolerance, alpha, A->root, B->root, &(C->root), timer, sgemm);
}

/** @brief Multiply a node by a scalar.
 *
 * The node and its children nodes are multiplied by the scalar.
 *
 * @param beta The scalar.
 * @param node The node.
 */
void
spamm_naive_multiply_scalar (const float beta, struct spamm_naive_node_t *node)
{
  int i;

  if (node == NULL) { return; }

  if (node->data != NULL)
  {
    for (i = 0; i < node->blocksize*node->blocksize; i++)
    {
      node->data[i] *= beta;
    }
  }

  else
  {
    for (i = 0; i < 4; i++)
    {
      spamm_naive_multiply_scalar(beta, node->child[i]);
    }
  }
}

/** Multiply to matrices, i.e. \f$ C = \alpha A \times B + \beta C\f$.
 *
 * @param tolerance The SpAMM tolerance of this product.
 * @param alpha The paramater \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param beta The paramater \f$\beta\f$.
 * @param C The matrix \f$C\f$.
 * @param timer The timer to use.
 * @param sgemm The external sgemm function to use.
 */
void
spamm_naive_multiply_matrix (const float tolerance,
    const float alpha,
    struct spamm_naive_node_t *node_A,
    struct spamm_naive_node_t *node_B,
    struct spamm_naive_node_t **node_C,
    struct spamm_timer_t *timer,
    void (*sgemm) (char *, char *, int *, int *, int *, float *, float *, int *, float *, int *, float *, float *, int *))
{
  float beta = 1.0;
  int i, j, k;

  if (node_A == NULL || node_B == NULL) { return; }

  /* We have to allocate a new C block a tier up. */
  if (*node_C == NULL)
  {
    printf("[%s:%i] node_C should not be NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  /* Multiply matrix blocks. */
  if (node_A->data != NULL && node_B->data != NULL)
  {
    if (node_A->norm*node_B->norm > tolerance)
    {
      if ((*node_C)->data == NULL)
      {
        (*node_C)->data = calloc((*node_C)->blocksize*(*node_C)->blocksize, sizeof(float));
      }
      sgemm("N", "N", &(node_A->blocksize), &(node_A->blocksize), &(node_A)->blocksize,
          (float*) &alpha, node_A->data, &node_A->blocksize, node_B->data,
          &node_A->blocksize, &beta, (*node_C)->data, &node_A->blocksize);
    }
  }

  else
  {
    /* Recurse. */
    for (i = 0; i < 2; i++) {
      for (j = 0; j < 2; j++) {
        for (k = 0; k < 2; k++)
        {
          if (node_A->child[spamm_index_row_major(i, k, 2, 2)] != NULL && node_B->child[spamm_index_row_major(k, j, 2, 2)] != NULL)
          {
            if (node_A->child[spamm_index_row_major(i, k, 2, 2)]->norm *
                node_B->child[spamm_index_row_major(k, j, 2, 2)]->norm > tolerance)
            {
              /* Create a new C node if necessary. */
              if ((*node_C)->child[spamm_index_row_major(i, j, 2, 2)] == NULL)
              {
                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)] = spamm_naive_new_node((*node_C)->tier+1);
                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)]->M_lower = (*node_C)->M_lower+((*node_C)->M_upper-(*node_C)->M_lower)/2*i;
                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)]->M_upper = (*node_C)->M_lower+((*node_C)->M_upper-(*node_C)->M_lower)/2*(i+1);
                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)]->N_lower = (*node_C)->N_lower+((*node_C)->N_upper-(*node_C)->N_lower)/2*j;
                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)]->N_upper = (*node_C)->N_lower+((*node_C)->N_upper-(*node_C)->N_lower)/2*(j+1);
              }

              spamm_naive_multiply_matrix(tolerance,
                  alpha,
                  node_A->child[spamm_index_row_major(i, k, 2, 2)],
                  node_B->child[spamm_index_row_major(k, j, 2, 2)],
                  &((*node_C)->child[spamm_index_row_major(i, j, 2, 2)]),
                  timer,
                  sgemm);
            }
          }
        }
      }
    }
  }
}

/** @brief Print a naive SpAMM matrix.
 *
 * @param A The matrix.
 */
void
spamm_naive_print (const struct spamm_naive_t *A)
{
  int i, j;

  if (A == NULL)
  {
    printf("A = (null)\n");
  }

  else
  {
    for (i = 0; i < A->M; i++) {
      for (j = 0; j < A->N; j++)
      {
        printf(" % 1.2e", spamm_naive_get(i, j, A));
      }
      printf("\n");
    }
  }
}
