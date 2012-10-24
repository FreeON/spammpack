/** @file */

#ifndef __SPAMM_CHUNK_H
#define __SPAMM_CHUNK_H

#include "spamm_types.h"

#include <stdlib.h>

spamm_chunk_t *
spamm_new_chunk (const unsigned int number_dimensions,
    const unsigned int N_contiguous);

size_t
spamm_chunk_align (const size_t address,
    const size_t alignment);

unsigned int *
spamm_chunk_get_number_dimensions (spamm_chunk_t *chunk);

unsigned int *
spamm_chunk_get_N_contiguous (spamm_chunk_t *chunk);

unsigned int *
spamm_chunk_get_N_lower (spamm_chunk_t *chunk);

unsigned int *
spamm_chunk_get_N_upper (spamm_chunk_t *chunk);

float *
spamm_chunk_get_matrix (spamm_chunk_t *chunk);

float *
spamm_chunk_get_matrix_dilated (spamm_chunk_t *chunk);

float *
spamm_chunk_get_norm (spamm_chunk_t *chunk);

float *
spamm_chunk_get_norm2 (spamm_chunk_t *chunk);

unsigned int
spamm_chunk_get_number_tiers (spamm_chunk_t *chunk);

size_t
spamm_chunk_get_size (const unsigned int number_dimensions,
    const unsigned int N_contiguous,
    unsigned int **number_dimension_pointer,
    unsigned int **N_contiguous_pointer,
    unsigned int **N_lower_pointer,
    unsigned int **N_upper_pointer,
    float **A_pointer,
    float **A_dilated_pointer,
    float **norm_pointer,
    float **norm2_pointer);

void
spamm_chunk_print (spamm_chunk_t *chunk);

void
spamm_delete_chunk (spamm_chunk_t **chunk);

#endif
