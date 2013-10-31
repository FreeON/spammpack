/** @file */

#ifndef __SPAMM_CHUNK_H
#define __SPAMM_CHUNK_H

#include "spamm_types.h"

#include <stdlib.h>

spamm_norm_t
spamm_chunk_multiply_scalar (const float alpha,
    spamm_chunk_t *chunk,
    double *const flop,
    double *const memop);

spamm_norm_t
spamm_chunk_multiply (const spamm_norm_t tolerance,
    const float alpha,
    const spamm_chunk_t *const chunk_A,
    const spamm_chunk_t *const chunk_B,
    spamm_chunk_t *const chunk_C,
    sgemm_func sgemm,
    double *const flop,
    double *const memop);

size_t
spamm_chunk_pad (const size_t address,
    const size_t alignment);

unsigned int *
spamm_chunk_get_number_dimensions (const spamm_chunk_t *const chunk);

unsigned int *
spamm_chunk_get_number_tiers (const spamm_chunk_t *const chunk);

unsigned int *
spamm_chunk_get_use_linear_tree (const spamm_chunk_t *const chunk);

unsigned int *
spamm_chunk_get_N (const spamm_chunk_t *const chunk);

unsigned int *
spamm_chunk_get_N_lower (const spamm_chunk_t *const chunk);

unsigned int *
spamm_chunk_get_N_upper (const spamm_chunk_t *const chunk);

unsigned int
spamm_chunk_get_N_contiguous (const spamm_chunk_t *const chunk);

float *
spamm_chunk_get_matrix (const spamm_chunk_t *const chunk);

float *
spamm_chunk_get_matrix_dilated (const spamm_chunk_t *const chunk);

unsigned int
spamm_chunk_get_number_norm_entries (const spamm_chunk_t *const chunk);

spamm_norm_t *
spamm_chunk_get_norm (const spamm_chunk_t *const chunk);

spamm_norm_t *
spamm_chunk_get_tier_norm (const unsigned int tier,
    const spamm_chunk_t *const chunk);

spamm_norm_t *
spamm_chunk_get_norm2 (const spamm_chunk_t *const chunk);

spamm_norm_t *
spamm_chunk_get_tier_norm2 (const unsigned int tier,
    const spamm_chunk_t *const chunk);

unsigned int
spamm_chunk_matrix_index (const unsigned int number_dimensions,
    const short use_linear_tree,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    const unsigned int *const i);

unsigned int
spamm_chunk_norm_index (const unsigned int tier,
    const unsigned int *const i,
    const spamm_chunk_t *const chunk);

size_t
spamm_chunk_get_size (const unsigned int number_dimensions,
    const short use_linear_tree,
    unsigned int *number_tiers,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    unsigned int **N_pointer,
    unsigned int **N_lower_pointer,
    unsigned int **N_upper_pointer,
    float **A_pointer,
    float **A_dilated_pointer,
    spamm_norm_t **norm_pointer,
    spamm_norm_t **norm2_pointer);

spamm_norm_t
spamm_chunk_fix (spamm_chunk_t *const chunk,
    double *const flop,
    double *const memop);

#endif
