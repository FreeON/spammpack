/** @file
 *
 * header file for spamm_introspection.c
 */

#ifndef __SPAMM_INTROSPECTION_H
#define __SPAMM_INTROSPECTION_H

#include "spamm_types.h"

unsigned int
spamm_get_number_dimensions (const struct spamm_matrix_t *const A);

unsigned int *
spamm_get_N (const struct spamm_matrix_t *const A);

#endif
