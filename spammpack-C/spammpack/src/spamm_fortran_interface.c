/** @file */

#include "spamm.h"

#include <stdio.h>
#include <stdlib.h>

void
spamm_new_chunk_interface (unsigned int *number_dimensions,
    unsigned int *N_contiguous,
    void **chunk)
{
  *chunk = spamm_new_chunk(*number_dimensions, *N_contiguous);
}

void
spamm_chunk_get_N_lower_interface (uint32_t **N_lower, spamm_chunk_t **chunk)
{
  *N_lower = spamm_chunk_get_N_lower(*chunk);
}

void
spamm_chunk_get_N_upper_interface (uint32_t **N_upper, spamm_chunk_t **chunk)
{
  *N_upper = spamm_chunk_get_N_upper(*chunk);
}

void
spamm_chunk_print_interface (spamm_chunk_t **chunk)
{
  spamm_chunk_print(*chunk);
}
