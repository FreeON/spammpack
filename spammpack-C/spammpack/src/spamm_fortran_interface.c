/** @file */

#include "config.h"
#include "spamm.h"

/** Fortran interface for spamm_trace().
 */
void
spamm_trace_interface (float *const trace,
    const struct spamm_matrix_t **const A)
{
  *trace = spamm_trace(*A);
}

/** Fortran interface for spamm_multiply().
 */
void
spamm_multiply_interface (const spamm_norm_t *const tolerance,
    const float *const alpha,
    const struct spamm_matrix_t **const A,
    const struct spamm_matrix_t **const B,
    const float *const beta,
    struct spamm_matrix_t **const C,
    double *const flop,
    double *const mop)
{
  spamm_multiply(*tolerance, *alpha, *A, *B, *beta, *C, SGEMM, flop, mop);
}

/** Fortran interface for spamm_multiply_scalar().
 */
void
spamm_multiply_scalar_interface (const float *const alpha,
    struct spamm_matrix_t **const A,
    double *const flop,
    double *const mop)
{
  spamm_multiply_scalar(*alpha, *A, flop, mop);
}

/** Fortran interface for spamm_get().
 */
void
spamm_get_interface (const unsigned int *const i,
    const struct spamm_matrix_t **const A,
    float *const Aij)
{
  *Aij = spamm_get(i, *A);
}

/** Fortran interface for spamm_chunk_multiply_scalar().
 */
void
spamm_chunk_multiply_scalar_interface (const float *const alpha,
    spamm_chunk_t **chunk, spamm_norm_t *const norm2,
    double *const flop,
    double *const mop)
{
  *norm2 = spamm_chunk_multiply_scalar(*alpha, *chunk, flop, mop);
}

/** Fortran interface for spamm_chunk_multiply().
 */
void
spamm_chunk_multiply_interface (const spamm_norm_t *const tolerance,
    const float *const alpha,
    spamm_chunk_t **chunk_A,
    spamm_chunk_t **chunk_B,
    spamm_chunk_t **chunk_C,
    spamm_norm_t *const norm2,
    double *const flop,
    double *const mop)
{
  *norm2 = spamm_chunk_multiply(*tolerance, *alpha, *chunk_A, *chunk_B,
      *chunk_C, SGEMM, flop, mop);
}

/** Fortran interface for spamm_convert_dense_to_spamm().
 */
void
spamm_convert_dense_to_spamm_interface (const unsigned int *const number_dimensions,
    const unsigned int *const N,
    const unsigned int *const chunk_tier,
    const unsigned int *const use_linear_tree,
    const float *const A_dense,
    struct spamm_matrix_t **const A)
{
  *A = spamm_convert_dense_to_spamm(*number_dimensions, N, *chunk_tier,
      *use_linear_tree, column_major, A_dense);
}

/** Fortran interface for spamm_convert_spamm_to_dense().
 */
void
spamm_convert_spamm_to_dense_interface (float **const A_dense,
    struct spamm_matrix_t **const A)
{
  *A_dense = spamm_convert_spamm_to_dense(*A);
}

/** Fortran interface for spamm_new_chunk().
 *
 * @param number_dimensions The number of dimensions.
 * @param N_contiguous The size of the contiguous submatrix.
 * @param chunk[out] The chunk.
 */
void
spamm_new_chunk_interface (unsigned int *number_dimensions,
    unsigned int *use_linear_tree,
    unsigned int *N,
    unsigned int *N_lower,
    unsigned int *N_upper,
    spamm_chunk_t **chunk)
{
  *chunk = spamm_new_chunk(*number_dimensions, *use_linear_tree, N, N_lower, N_upper);
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
spamm_chunk_get_norm_interface (spamm_norm_t **norm, spamm_chunk_t **chunk)
{
  *norm = spamm_chunk_get_norm(*chunk);
}

/** Fortran interface wrapper for spamm_chunk_get_norm2().
 *
 * @param norm2 The norm2 vector.
 * @param chunk The chunk.
 */
void
spamm_chunk_get_norm2_interface (spamm_norm_t **norm2, spamm_chunk_t **chunk)
{
  *norm2 = spamm_chunk_get_norm2(*chunk);
}

void
spamm_print_chunk_interface (spamm_chunk_t **chunk)
{
  spamm_print_chunk(*chunk);
}

/** Fortran interface for spamm_new().
 */
void
spamm_new_interface (const unsigned int *const number_dimensions,
    const unsigned int *const N,
    const unsigned int *const chunk_tier,
    const unsigned int *const use_linear_tree,
    struct spamm_matrix_t **const A)
{
  *A = spamm_new(*number_dimensions, N, *chunk_tier, *use_linear_tree);
}

/** Fortran interface for spamm_copy().
 */
void
spamm_copy_interface (struct spamm_matrix_t **const A,
    const float *alpha,
    const struct spamm_matrix_t **const B,
    double *const flop,
    double *const mop)
{
  spamm_copy(A, *alpha, *B, flop, mop);
}

/** Fortran interface for spamm_add(), @f$ C \leftarrow \alpha A + \beta B
 * @f$.
 */
void
spamm_add_interface (const float *const alpha,
    const struct spamm_matrix_t **const A,
    const float *const beta,
    const struct spamm_matrix_t **const B,
    struct spamm_matrix_t **const C,
    double *const flop,
    double *const mop)
{
  spamm_copy(C, 1.0, *A, flop, mop);
  spamm_add(*alpha, *C, *beta, *B, flop, mop);
}
