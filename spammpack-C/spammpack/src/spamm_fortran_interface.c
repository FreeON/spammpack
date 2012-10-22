/** @file */

#include "spamm.h"

void
spamm_new_chunk_interface (unsigned int *number_dimensions,
    unsigned int *N_contiguous,
    void **chunk)
{
  *chunk = spamm_new_chunk(*number_dimensions, *N_contiguous);
}

void
spamm_chunk_get_N_lower_interface (unsigned int **N_lower, spamm_chunk_t **chunk)
{
  *N_lower = spamm_chunk_get_N_lower(*chunk);
}

void
spamm_chunk_get_N_upper_interface (unsigned int **N_upper, spamm_chunk_t **chunk)
{
  *N_upper = spamm_chunk_get_N_upper(*chunk);
}

/** Fortran interface wrapper for spamm_chunk_get_matrix().
 *
 * @param A The matrix.
 * @param chunk The chunk.
 */
void
spamm_chunk_get_matrix_interface (float **A, spamm_chunk_t **chunk)
{
  *A = spamm_chunk_get_matrix(*chunk);
}

/** Fortran interface wrapper for spamm_chunk_get_matrix_dilated().
 *
 * @param A The matrix.
 * @param chunk The chunk.
 */
void
spamm_chunk_get_matrix_dilated_interface (float **A, spamm_chunk_t **chunk)
{
  *A = spamm_chunk_get_matrix_dilated(*chunk);
}

void
spamm_chunk_print_interface (spamm_chunk_t **chunk)
{
  spamm_chunk_print(*chunk);
}
