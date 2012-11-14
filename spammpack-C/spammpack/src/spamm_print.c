#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/** Print a dense matrix.
 *
 * @param M The number of rows.
 * @param N The number of columns.
 * @param type The storage type of the dense matrix.
 * @param A The dense matrix.
 */
void
spamm_print_dense (const unsigned int M, const unsigned int N,
    const enum spamm_layout_t type, const float *A)
{
  unsigned int i, j;

  switch (type)
  {
    case column_major:
      for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++)
        {
          printf(" % 1.2e", A[spamm_index_column_major(i, j, M, N)]);
        }
        printf("\n");
      }
      break;

    case row_major:
      for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++)
        {
          printf(" % 1.2e", A[spamm_index_row_major(i, j, M, N)]);
        }
        printf("\n");
      }
      break;

    default:
      printf("error\n");
      exit(1);
      break;
  }
}

/** Print a SpAMM chunk.
 *
 * @param chunk The chunk.
 */
void
spamm_print_chunk (spamm_chunk_t *const chunk)
{
  unsigned int i, j;
  unsigned int N_contiguous;
  float *matrix;

  printf("chunk (%p): ", chunk);
  printf("ndim = %u", *spamm_chunk_get_number_dimensions(chunk));
  printf(", ntiers = %u", *spamm_chunk_get_number_tiers(chunk));
  printf(", lintree = %u", *spamm_chunk_get_use_linear_tree(chunk));

  N_contiguous = spamm_chunk_get_N_contiguous(chunk);
  matrix = spamm_chunk_get_matrix(chunk);
  printf(", matrix =\n");
  spamm_print_dense(N_contiguous, N_contiguous, column_major, matrix);
}

/** Print a recursive tree node and all nodes underneath.
 *
 * @param node The tree node.
 */
void
spamm_recursive_print_node (const struct spamm_recursive_node_t *const node,
    const unsigned int number_dimensions,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const short use_linear_tree)
{
  unsigned int i;
  int dim;

  unsigned int *new_N_lower;
  unsigned int *new_N_upper;

  if (node == NULL) { return; }

  printf("node (tier = %u, %p): N = {", tier, node);
  for (dim = 0; dim < number_dimensions; dim++)
  {
    printf(" [ %u --> %u ]", N_lower[dim], N_upper[dim]);
  }
  printf(" }, norm = %1.2e, ", node->norm);
  if (tier == chunk_tier)
  {
    spamm_print_chunk(node->tree.chunk);
  }

  else if (tier == chunk_tier)
  {
    return;
  }

  else
  {
    if (node->tree.child == NULL)
    {
      printf("no children");
    }

    else
    {
      printf("{");
      for (i = 0; i < number_dimensions*number_dimensions; i++)
      {
        printf(" %p", node->tree.child[i]);
      }
      printf(" }");
    }
    printf("\n");

    for (i = 0; i < ipow(2, number_dimensions); i++)
    {
      new_N_lower = malloc(number_dimensions*sizeof(unsigned int));
      new_N_upper = malloc(number_dimensions*sizeof(unsigned int));

      for (dim = 0; dim < number_dimensions; dim++)
      {
        new_N_lower[dim] = N_lower[dim]+(N_upper[dim]-N_lower[dim])/2*(i & (1 << dim) ? 1 : 0);
        new_N_upper[dim] = new_N_lower[dim]+(N_upper[dim]-N_lower[dim])/2;
      }

      spamm_recursive_print_node(node->tree.child[i], number_dimensions,
          new_N_lower, new_N_upper, tier+1, chunk_tier, use_linear_tree);

      free(new_N_lower);
      free(new_N_upper);
    }
  }
}

/** Print a SpAMM matrix, printing the tree structure.
 *
 * @param A The matrix.
 */
void
spamm_print_tree (const struct spamm_matrix_t *const A)
{
  int dim;
  unsigned int tier;

  unsigned int *N_lower;
  unsigned int *N_upper;

  assert(A != NULL);

  printf("root node: ");
  printf("ndim = %u, ", A->number_dimensions);
  printf("N = {");
  for (dim = 0; dim < A->number_dimensions; dim++)
  {
    printf(" %u", A->N[dim]);
  }
  printf(" }, ");
  printf("N_padded = %u, ", A->N_padded);
  printf("depth = %u, ", A->depth);
  printf("chunk_tier = %u, ", A->chunk_tier);
  printf("use_linear_tree = %u, ", A->use_linear_tree);
  printf("kernel_tier = %u, ", A->kernel_tier);
  printf("chunk_tier = %u\n", A->chunk_tier);

  if (A->chunk_tier == 0)
  {
    spamm_print_chunk(A->tree.chunk);
  }
  else
  {
    N_lower = calloc(A->number_dimensions, sizeof(unsigned int));
    N_upper = calloc(A->number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < A->number_dimensions; dim++)
    {
      N_upper[dim] = A->N_padded;
    }

    spamm_recursive_print_node(A->tree.recursive_tree, A->number_dimensions,
        N_lower, N_upper, 0, A->chunk_tier, A->use_linear_tree);

    free(N_lower);
    free(N_upper);
  }
}

/** Print a SpAMM matrix in the same format a dense matrix is printed.
 *
 * @param A The matrix.
 */
void
spamm_print (const struct spamm_matrix_t *A)
{
  unsigned int *i;

  assert(A != NULL);

  i = calloc(A->number_dimensions, sizeof(unsigned int));
  switch(A->number_dimensions)
  {
    case 2:
      for (i[0] = 0; i[0] < A->N[0]; i[0]++) {
        for (i[1] = 0; i[1] < A->N[1]; i[1]++)
        {
          printf(" % 1.2e", spamm_get(i, A));
        }
        printf("\n");
      }
      break;

    default:
      SPAMM_FATAL("not implemented\n");
  }
}
