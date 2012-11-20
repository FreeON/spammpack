#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <stdio.h>
#include <stdint.h>
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
  int dim;
  unsigned int i, j;
  unsigned int number_dimensions;
  unsigned int N_contiguous;
  unsigned int use_linear_tree;
  unsigned int number_tiers;

  unsigned int *N;
  unsigned int *N_lower;
  unsigned int *N_upper;

  float *A;
  float *norm;
  float *norm2;

  if (chunk == NULL) { return; }

  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
  number_tiers = *spamm_chunk_get_number_tiers(chunk);
  use_linear_tree = *spamm_chunk_get_use_linear_tree(chunk);
  N_contiguous = spamm_chunk_get_N_contiguous(chunk);

  printf("chunk (%p): ", chunk);
  printf("ndim = %u", number_dimensions);
  printf(", ntiers = %u", number_tiers);
  printf(", N_cont = %u", N_contiguous);
  printf(", lintree = %u", use_linear_tree);
  printf("\n");
  N = spamm_chunk_get_N(chunk);
  printf("N = [");
  for (dim = 0; dim < number_dimensions; dim++)
  {
    printf(" %u", N[dim]);
    if (dim+1 < number_dimensions)
    {
      printf(",");
    }
  }
  printf(" ]\n");
  N_lower = spamm_chunk_get_N_lower(chunk);
  N_upper = spamm_chunk_get_N_upper(chunk);
  for (dim = 0; dim < number_dimensions; dim++)
  {
    printf("N[%u] = [ %u, %u ]\n", dim, N_lower[dim], N_upper[dim]);
  }

  for (i = 0; i < number_tiers; i++)
  {
    printf("&norm[%u] = %p (%p)\n", i, (void*) ((intptr_t) spamm_chunk_get_tier_norm(i, chunk) - (intptr_t) chunk),
        spamm_chunk_get_tier_norm(i, chunk));
  }

  for (i = 0; i < number_tiers; i++)
  {
    printf("&norm2[%u] = %p (%p)\n", i, (void*) ((intptr_t) spamm_chunk_get_tier_norm2(i, chunk) - (intptr_t) chunk),
        spamm_chunk_get_tier_norm2(i, chunk));
  }

  printf("&matrix = %p (%p)\n", (void*) ((intptr_t) spamm_chunk_get_matrix(chunk) - (intptr_t) chunk),
      spamm_chunk_get_matrix(chunk));

  for (i = 0; i < number_tiers; i++)
  {
    norm = spamm_chunk_get_tier_norm(i, chunk);
    printf("norm[%u] =", i);
    for (j = 0; j < ipow(4, i); j++)
    {
      printf(" %1.2f", norm[j]);
    }
    printf("\n");
  }

  for (i = 0; i < number_tiers; i++)
  {
    norm2 = spamm_chunk_get_tier_norm2(i, chunk);
    printf("norm2[%u] =", i);
    for (j = 0; j < ipow(4, i); j++)
    {
      printf(" %1.2f", norm2[j]);
    }
    printf("\n");
  }

  printf("matrix:\n");
  A = spamm_chunk_get_matrix(chunk);
  spamm_print_dense(N_contiguous, N_contiguous, column_major, A);
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
