/** @file
 *
 * Header file for generate_matrix.c
 */

#ifndef __GENERATE_MATRIX_H
#define __GENERATE_MATRIX_H

float *
generate_matrix (const unsigned int number_dimensions,
    const short is_sparse,
    unsigned int **const N);

struct spamm_matrix_t *
create_spamm_from_dense (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int chunk_tier,
    const short use_linear_tree,
    const float *const A_dense);

int
compare_spamm_to_dense (const struct spamm_matrix_t *const A,
    const float *const A_dense);

#endif
