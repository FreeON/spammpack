/** @file
 *
 * Header file for generate_matrix.c
 */

#ifndef __TEST_H
#define __TEST_H

#include <spamm_types.h>

unsigned int *
generate_shape (const unsigned int number_dimensions,
    const short is_square);

float *
generate_matrix (const unsigned int number_dimensions,
    const short is_sparse,
    const unsigned int *const N);

int
compare_spamm_to_dense (const struct spamm_matrix_t *const A,
    const float *const A_dense,
    const double abs_tolerance);

#endif
