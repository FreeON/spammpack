#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/** Return the number of non-zero elements in a matrix.
 *
 * @param A The matrix.
 *
 * @return The number of non-zero elements.
 */
unsigned int
spamm_number_nonzero (const struct spamm_matrix_t *const A)
{
  unsigned int *i;
  unsigned int result = 0;

  i = calloc(A->number_dimensions, sizeof(unsigned int));

  switch(A->number_dimensions)
  {
    case 1:
      for(i[0] = 0; i[0] < A->N[0]; i[0]++)
      {
        if(spamm_get(i, A) != 0.0)
        {
          result++;
        }
      }
      break;

    case 2:
      for(i[0] = 0; i[0] < A->N[0]; i[0]++) {
        for(i[1] = 0; i[1] < A->N[1]; i[1]++)
        {
          if(spamm_get(i, A) != 0.0)
          {
            result++;
          }
        }
      }
      break;

    case 3:
      for(i[0] = 0; i[0] < A->N[0]; i[0]++) {
        for(i[1] = 0; i[1] < A->N[1]; i[1]++) {
          for(i[2] = 0; i[2] < A->N[2]; i[2]++)
          {
            if(spamm_get(i, A) != 0.0)
            {
              result++;
            }
          }
        }
      }
      break;

    default:
      SPAMM_FATAL("not implemented\n");
  }

  free(i);

  return result;
}

/** Count the number of nodes in the matrix.
 *
 * @param tier The current tier.
 * @param chunk_tier The tier at which SpAMM chunks are stored.
 * @param number_dimensions The number of dimensions.
 * @param number_recursive_nodes [out] The number of recursive nodes.
 * @param number_chunks [out] The number of chunks.
 * @param node The node.
 */
void
spamm_recursive_number_nodes (const unsigned int tier,
    const unsigned int chunk_tier,
    const unsigned int number_dimensions,
    unsigned int *const number_recursive_nodes,
    unsigned int *const number_chunks,
    const struct spamm_recursive_node_t *const node)
{
  unsigned int i;

  if(node == NULL) { return; }

  (*number_recursive_nodes)++;

  if(tier == chunk_tier)
  {
    (*number_chunks)++;
  }

  else
  {
    for(i = 0; i < ipow(2, number_dimensions); i++)
    {
      spamm_recursive_number_nodes(tier+1, chunk_tier, number_dimensions,
          number_recursive_nodes, number_chunks, node->tree.child[i]);
    }
  }
}

/** Count the number of nodes in the matrix.
 *
 * @param number_recursive_nodes [out] The number of recursive nodes.
 * @param number_chunks [out] The number of chunks.
 * @param [in] A The matrix.
 */
void
spamm_number_nodes (unsigned int *const number_recursive_nodes,
    unsigned int *const number_chunks,
    const struct spamm_matrix_t *const A)
{
  unsigned int i;

  *number_recursive_nodes = 0;
  *number_chunks = 0;

  if(A == NULL) { return; }

  (*number_recursive_nodes)++;

  if(A->chunk_tier == 0)
  {
    (*number_chunks)++;
  }

  else
  {
    for(i = 0; i < ipow(2, A->number_dimensions); i++)
    {
      spamm_recursive_number_nodes(1, A->chunk_tier, A->number_dimensions,
          number_recursive_nodes, number_chunks, A->tree.recursive_tree);
    }
  }
}

/** Print out some information on the matrix.
 *
 * @param A The matrix.
 */
void
spamm_print_info (const struct spamm_matrix_t *const A)
{
  int dim;
  unsigned int i;
  unsigned int N_matrix;
  unsigned int N_contiguous;
  unsigned int number_recursive_nodes;
  unsigned int number_recursive_nodes_full;
  unsigned int number_chunks;
  unsigned int number_chunks_full;

  assert(A != NULL);

  for(dim = 0, N_matrix = 1; dim < A->number_dimensions; dim++)
  {
    N_matrix *= A->N[dim];
  }

  for(i = 0, N_contiguous = A->N_padded; i < A->chunk_tier; i++)
  {
    N_contiguous >>= 1;
  }

  spamm_number_nodes(&number_recursive_nodes, &number_chunks, A);

  for(i = 0, number_recursive_nodes_full = 0; i <= A->chunk_tier; i++)
  {
    number_recursive_nodes_full += ipow(ipow(2, A->number_dimensions), i);
  }
  number_chunks_full = ipow(ipow(2, A->number_dimensions), A->chunk_tier);

  printf("number_dimensions = %u", A->number_dimensions);
  printf(", recursive nodes = %u (out of %u)", number_recursive_nodes, number_recursive_nodes_full);
  printf(", chunks = %u (out of %u)", number_chunks, number_chunks_full);
  printf(", N = {");
  for(dim = 0; dim < A->number_dimensions; dim++)
  {
    printf(" %u", A->N[dim]);
    if(dim+1 < A->number_dimensions)
    {
      printf(",");
    }
  }
  printf(" }");
  printf(", N_padded = %u", A->N_padded);
  printf(", chunk_tier = %u", A->chunk_tier);
  printf(", N_contiguous = %u", N_contiguous);
  printf(", use_linear_tree = %u", A->use_linear_tree);
  printf(", nnzero = %u (%1.2f%% sparsity)", spamm_number_nonzero(A), 100*(1-(float) spamm_number_nonzero(A)/(float) N_matrix));
  printf("\n");
}
