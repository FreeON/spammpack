/** @file */

#include "spamm.h"

/** Fortran interface for spamm_new_chunk().
 *
 * @param number_dimensions The number of dimensions.
 * @param N_contiguous The size of the contiguous submatrix.
 * @param chunk[out] The chunk.
 */
void
spamm_new_chunk_interface (unsigned int *number_dimensions,
    unsigned int *N_block,
    unsigned int *N,
    unsigned int *N_lower,
    unsigned int *N_upper,
    spamm_chunk_t **chunk)
{
  *chunk = spamm_new_chunk(*number_dimensions, *N_block, N, N_lower, N_upper);
}

void
spamm_delete_chunk_interface (spamm_chunk_t **chunk)
{
  spamm_delete_chunk(chunk);
}

void
spamm_chunk_get_number_dimensions_interface (unsigned int **number_dimensions, spamm_chunk_t **chunk)
{
  *number_dimensions = spamm_chunk_get_number_dimensions(*chunk);
}

void
spamm_chunk_get_N_lower_interface (unsigned int **N_lower, spamm_chunk_t **chunk)
{
  *N_lower = spamm_chunk_get_N_lower(*chunk);
}

void
spamm_chunk_get_N_contiguous_interface (unsigned int *N_contiguous, spamm_chunk_t **chunk)
{
  *N_contiguous = spamm_chunk_get_N_contiguous(*chunk);
}

void
spamm_chunk_get_number_tiers_interface (unsigned int **number_tiers, spamm_chunk_t **chunk)
{
  *number_tiers = spamm_chunk_get_number_tiers(*chunk);
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

/** Fortran interface wrapper for spamm_chunk_get_norm().
 *
 * @param norm The norm vector.
 * @param chunk The chunk.
 */
void
spamm_chunk_get_norm_interface (float **norm, spamm_chunk_t **chunk)
{
  *norm = spamm_chunk_get_norm(*chunk);
}

/** Fortran interface wrapper for spamm_chunk_get_norm2().
 *
 * @param norm2 The norm2 vector.
 * @param chunk The chunk.
 */
void
spamm_chunk_get_norm2_interface (float **norm2, spamm_chunk_t **chunk)
{
  *norm2 = spamm_chunk_get_norm2(*chunk);
}

void
spamm_chunk_print_interface (spamm_chunk_t **chunk)
{
  spamm_chunk_print(*chunk);
}
