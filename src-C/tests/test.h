/** @file
 *
 * Header file for generate_matrix.c
 */

#ifndef __TEST_H
#define __TEST_H

#include <spamm_types.h>

/** The number of matrix types. */
#define NUMBER_MATRIX_TYPES 5

/** The matrix types for the tests. */
enum matrix_t
{
  /** A full matrix. */
  full,

  /** A sparse matrix that is diagonally banded. */
  diagonally_banded,

  /** A diagonally dominant matrix with exponential decay away from the
   * diagonal. */
  exponential_decay,

  /** A sparse matrix that has few randomly placed non-zero elements. */
  sparse_random,

  /** A matrix with just a handfull of non-zero elements. */
  hyper_sparse
};

const char *const
get_matrix_type_name (const enum matrix_t matrix_type);

char *
print_matrix_types ();

enum matrix_t
parse_matrix_type (const char *const type_name);

unsigned int *
generate_shape (const unsigned int number_dimensions,
    const short is_square,
    const unsigned int N_mean,
    const unsigned int N_width);

float *
generate_matrix_float (const unsigned int number_dimensions,
    const enum matrix_t matrix_type,
    const unsigned int *const N,
    const float gamma);

double *
generate_matrix_double (const unsigned int number_dimensions,
    const enum matrix_t matrix_type,
    const unsigned int *const N,
    const double gamma);

int
compare_spamm_to_dense_float (const struct spamm_matrix_t *const A,
    const float *const A_dense,
    const double abs_tolerance);

int
compare_spamm_to_dense_double (const struct spamm_matrix_t *const A,
    const double *const A_dense,
    const double abs_tolerance);

#endif
