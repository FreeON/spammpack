/** @file
 *
 * Get information on a matrix or a node.
 */

#include "spamm_types.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <stdlib.h>

/** Get the number of dimensions of a matrix.
 *
 * @param A The matrix.
 *
 * @return The number of dimensions.
 */
unsigned int
spamm_get_number_dimensions (const struct spamm_matrix_t *const A)
{
  assert(A != NULL);

  return A->number_dimensions;
}

/** Get the shape of the matrix.
 *
 * @param A The matrix.
 *
 * @return An array of size number_dimension with the shape.
 */
unsigned int *
spamm_get_N (const struct spamm_matrix_t *const A)
{
  assert(A != NULL);

  return A->N;
}

/** Get the Frobenius norm of a matrix.
 *
 * @param A The matrix.
 *
 * @return The Frobenius norm.
 */
spamm_norm_t
spamm_get_norm (const struct spamm_matrix_t *const A)
{
  if(A == NULL) { return 0; }
  else if(A->recursive_tree == NULL) { return 0; }
  else
  {
    return A->recursive_tree->norm;
  }
}
