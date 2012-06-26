/** @file */

#include "spamm_list.h"

struct spamm_list_t *spamm_matrix_object;

/** Convert dense matrix to spamm matrix.
 *
 * @param A_dense The dense matrix.
 * @param B The SpAMM_Matrix object identifier.
 */
void
spamm_convert_dense_to_spamm_interface (const float *A_dense, int *B)
{
}

void
spamm_convert_dense_to_spamm_interface_ (const float *A_dense, int *B)
{
  spamm_convert_dense_to_spamm_interface(A_dense, B);
}

/** Create a new spamm matrix.
 *
 * @param A The matrix.
 */
void
spamm_new_interface (int *A)
{
}

void
spamm_new_interface_ (int *A)
{
  spamm_new_interface(A);
}
