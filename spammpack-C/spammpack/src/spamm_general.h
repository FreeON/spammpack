/** @file */

#ifndef __SPAMM_GENERAL_H
#define __SPAMM_GENERAL_H

#include "spamm_types.h"

unsigned int
ipow (unsigned int b, unsigned int n);

void *
spamm_allocate (const size_t size, const short zero_memory);

int
spamm_check (const struct spamm_matrix_t *A, const float rel_tolerance);

void
spamm_copy (struct spamm_matrix_t **A,
    const float beta,
    const struct spamm_matrix_t *const B);

void
spamm_chunk_copy (spamm_chunk_t **A,
    const float beta,
    spamm_chunk_t *B,
    const short use_linear_tree);

void
spamm_recursive_copy (struct spamm_recursive_node_t *const A,
    const float beta,
    const struct spamm_recursive_node_t *const B,
    const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const short use_linear_tree);

void
spamm_delete_chunk (spamm_chunk_t **chunk);

void
spamm_recursive_delete (const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int chunk_tier,
    struct spamm_recursive_node_t **node);

void
spamm_delete (struct spamm_matrix_t **A);

unsigned int
spamm_index_row_major (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N);

unsigned int
spamm_index_column_major (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N);

unsigned int
spamm_index_row_major_2 (const unsigned int number_dimensions,
    const unsigned int N,
    const unsigned int *const i);

unsigned int
spamm_index_column_major_2 (const unsigned int number_dimensions,
    const unsigned int N,
    const unsigned int *const i);

unsigned int
spamm_index_row_major_3 (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int *const i);

unsigned int
spamm_index_column_major_3 (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int *const i);

unsigned int
spamm_index_norm (const unsigned int i, const unsigned int j);

unsigned int
spamm_index_kernel_block (const unsigned int i, const unsigned int j, const enum spamm_layout_t layout);

unsigned int
spamm_index_kernel_block_transpose (const unsigned int i, const unsigned int j, const enum spamm_layout_t layout);

unsigned int
spamm_index_kernel_block_hierarchical (const unsigned int i_blocked,
    const unsigned int j_blocked, const unsigned int i_basic,
    const unsigned int j_basic, const enum spamm_layout_t layout);

unsigned int
spamm_index_kernel_block_transpose_hierarchical (const unsigned int i_block,
    const unsigned int j_block, const unsigned int i,
    const unsigned int j, const enum spamm_layout_t layout);

float
spamm_get (const unsigned int *const i, const struct spamm_matrix_t *A);

spamm_norm_t
spamm_get_norm (const struct spamm_matrix_t *const A);

unsigned int
spamm_index_linear (const unsigned int number_dimensions,
    const unsigned int *const i);

unsigned int
spamm_index_2D (const unsigned int i, const unsigned int j);

void
spamm_index_2D_to_ij (const unsigned int index, unsigned int *i, unsigned int *j);

unsigned int
spamm_index_3D_0kj (const unsigned int k, const unsigned int j);

unsigned int
spamm_index_3D_ik0 (const unsigned int i, const unsigned int k);

unsigned int
spamm_index_3D_i0j_to_2D (const unsigned int index_3D_i0j);

unsigned int
spamm_index_3D_ikj_to_k (const unsigned int index_3D_ikj);

double
spamm_linear_multiply (const float tolerance,
    const float alpha,
    const spamm_chunk_t *const chunk_A,
    const spamm_chunk_t *const chunk_B,
    spamm_chunk_t *const chunk_C,
    double *const flop);

void
spamm_recursive_multiply_scalar (const float alpha,
    struct spamm_recursive_node_t *A,
    const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const short use_linear_tree,
    double *const flop);

void
spamm_multiply_scalar (const float alpha,
    struct spamm_matrix_t *const A,
    double *const flop);

void
spamm_multiply (const float tolerance,
    const float alpha,
    const struct spamm_matrix_t *const A,
    const struct spamm_matrix_t *const B,
    const float beta,
    struct spamm_matrix_t *const C,
    sgemm_func sgemm,
    double *const flop);

void
spamm_add (const float alpha,
    struct spamm_matrix_t *const A,
    const float beta,
    const struct spamm_matrix_t *const B,
    double *const flop);

unsigned int
spamm_number_nonzero (const struct spamm_matrix_t *A);

spamm_chunk_t *
spamm_new_chunk (const unsigned int number_dimensions,
    const short use_linear_tree,
    const unsigned int *const N,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper);

struct spamm_matrix_t *
spamm_new (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int chunk_tier,
    const short use_linear_tree);

struct spamm_recursive_node_t *
spamm_recursive_new_node ();

void
spamm_print (const struct spamm_matrix_t *const A);

void
spamm_matlab_print (const struct spamm_matrix_t *const A);

void
spamm_print_info (const struct spamm_matrix_t *const A);

void
spamm_print_chunk (spamm_chunk_t *const chunk);

void
spamm_print_tree (const struct spamm_matrix_t *A);

void
spamm_print_dense (const unsigned int M, const unsigned int N,
    const enum spamm_layout_t type, const float *A);

void
spamm_set (const unsigned int *const i, const float Aij, struct spamm_matrix_t *A);

void
spamm_uint_to_bin_string (const unsigned int width, const unsigned int i, char *result);

char *
spamm_version ();

struct spamm_matrix_t *
spamm_convert_dense_to_spamm (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int chunk_tier,
    const short use_linear_tree,
    const enum spamm_layout_t dense_type,
    const float *const A_dense);

float *
spamm_convert_spamm_to_dense (const struct spamm_matrix_t *const A);

void
spamm_spectral_bounds (float *const a_min,
    float *const a_max,
    struct spamm_matrix_t *A);

void
spamm_sort_masked (const unsigned int length,
    unsigned int *const list,
    const unsigned int mask);

void spamm_sgemm (char * transA, char * transB,
    int *M, int *N, int *K,
    float *alpha, float *A, int *LDA, float *B, int *LDB,
    float *beta, float *C, int *LDC);

void spamm_dgemm (char * transA, char * transB,
    int *M, int *N, int *K,
    double *alpha, double *A, int *LDA, double *B, int *LDB,
    double *beta, double *C, int *LDC);

unsigned int
spamm_chunk_get_total_number_norms (const unsigned int number_tiers,
    const unsigned int number_dimensions);

float
spamm_trace (const struct spamm_matrix_t *const A);

void
spamm_omp_init ();

#endif
