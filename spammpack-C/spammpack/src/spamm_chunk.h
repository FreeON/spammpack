/** @file */

#ifndef __SPAMM_CHUNK_H
#define __SPAMM_CHUNK_H

#include "spamm_types.h"

#include <stdlib.h>

float
spamm_chunk_multiply_scalar (const float alpha,
    spamm_chunk_t *chunk);

float
spamm_chunk_multiply (const float tolerance,
    const float alpha,
    spamm_chunk_t *chunk_A,
    spamm_chunk_t *chunk_B,
    spamm_chunk_t *chunk_C,
    const unsigned int tier,
    const unsigned int contiguous_tier,
    const unsigned int depth,
    const unsigned int N_block,
    const unsigned int linear_index_A,
    const unsigned int linear_index_B,
    const unsigned int linear_index_C,
    sgemm_func sgemm);

void
spamm_chunk_add (const float alpha,
    spamm_chunk_t **A,
    const float beta,
    spamm_chunk_t *B);

void
spamm_chunk_copy (spamm_chunk_t **A,
    const float beta,
    spamm_chunk_t *B);

void
spamm_chunk_set (const unsigned int *const i,
    const float Aij,
    const unsigned int tier,
    const unsigned int contiguous_tier,
    const unsigned int depth,
    const unsigned int linear_index,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    spamm_chunk_t *chunk);

float
spamm_chunk_get (const unsigned int *i,
    const unsigned int tier,
    const unsigned int contiguous_tier,
    const unsigned int depth,
    const unsigned int linear_index,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    spamm_chunk_t *chunk);

spamm_chunk_t *
spamm_new_chunk (const unsigned int number_dimensions,
    const unsigned int N_block,
    const unsigned int *const N,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper);

size_t
spamm_chunk_pad (const size_t address,
    const size_t alignment);

unsigned int *
spamm_chunk_get_number_dimensions (spamm_chunk_t *chunk);

unsigned int *
spamm_chunk_get_N_block (spamm_chunk_t *chunk);

unsigned int *
spamm_chunk_get_N (spamm_chunk_t *chunk);

unsigned int *
spamm_chunk_get_N_lower (spamm_chunk_t *chunk);

unsigned int *
spamm_chunk_get_N_upper (spamm_chunk_t *chunk);

unsigned int
spamm_chunk_get_N_contiguous (spamm_chunk_t *chunk);

float *
spamm_chunk_get_matrix (spamm_chunk_t *chunk);

float *
spamm_chunk_get_matrix_dilated (spamm_chunk_t *chunk);

unsigned int
spamm_chunk_get_number_norm_entries (spamm_chunk_t *chunk);

float *
spamm_chunk_get_norm (spamm_chunk_t *chunk);

float *
spamm_chunk_get_tier_norm (const unsigned int tier,
    spamm_chunk_t *chunk);

float *
spamm_chunk_get_norm2 (spamm_chunk_t *chunk);

float *
spamm_chunk_get_tier_norm2 (const unsigned int tier,
    spamm_chunk_t *chunk);

unsigned int
spamm_chunk_get_number_tiers (spamm_chunk_t *chunk);

unsigned int
spamm_chunk_matrix_index (const unsigned int number_dimensions,
    const unsigned int N_block,
    const unsigned int *const N_lower,
    const unsigned int *const i);

size_t
spamm_chunk_get_size (const unsigned int number_dimensions,
    const unsigned int N_block,
    const unsigned int *const N,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    unsigned int **N_pointer,
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
