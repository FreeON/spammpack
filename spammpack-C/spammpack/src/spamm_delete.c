#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/** Delete a chunk. This function simply calls free().
 *
 * @param chunk The chunk to delete.
 */
void
spamm_delete_chunk (spamm_chunk_t **chunk)
{
  if(chunk == NULL) { return; }

  if(*chunk != NULL)
  {
    free(*chunk);
  }
  *chunk = NULL;
}

/** Delete a recursive matrix.
 *
 * @param A The recursive matrix root.
 */
void
spamm_recursive_delete (const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int chunk_tier,
    struct spamm_recursive_node_t **node)
{
  unsigned int i;

  if(*node == NULL) { return; }

  if(tier == chunk_tier)
  {
    spamm_delete_chunk(&(*node)->tree.chunk);
  }

  else
  {
    for(i = 0; i < ipow(2, number_dimensions); i++)
    {
      spamm_recursive_delete(number_dimensions, tier+1, chunk_tier, &(*node)->tree.child[i]);
    }
#ifdef _OPENMP
    omp_destroy_lock(&(*node)->lock);
#endif
    free((*node)->tree.child);
    (*node)->tree.child = NULL;
  }

  free(*node);
  *node = NULL;
}

/** Delete a matrix.
 *
 * @param A The matrix to delete.
 */
void
spamm_delete (struct spamm_matrix_t **A)
{
  assert(A != NULL);

  if(*A == NULL) { return; }

  spamm_recursive_delete((*A)->number_dimensions, 0, (*A)->chunk_tier, &(*A)->recursive_tree);

  /* Free memory. */
  free((*A)->N);
  free(*A);

  /* Reset pointer. */
  *A = NULL;
}
