/** @file */

#ifndef __SPAMM_GENERAL_H
#define __SPAMM_GENERAL_H

#include "spamm_types.h"

unsigned int
ipow (unsigned int b, unsigned int n);

void *
spamm_allocate (const size_t size, const short zero_memory);

int
spamm_check (const struct spamm_matrix_t *A, const float tolerance);

void
spamm_copy (struct spamm_matrix_t **A,
    const float beta,
    const struct spamm_matrix_t *const B);

void
spamm_recursive_copy (struct spamm_recursive_node_t **A,
    const float beta,
    const struct spamm_recursive_node_t *const B,
    const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int contiguous_tier,
    const short use_linear_tree);

void
spamm_hashed_copy (struct spamm_hashed_t **A,
    const float beta,
    const struct spamm_hashed_t *const B);

void
spamm_delete (struct spamm_matrix_t **A);

void
spamm_recursive_delete (const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int contiguous_tier,
    const short use_linear_tree,
    struct spamm_recursive_node_t **node);

void
spamm_delete_block (struct spamm_hashed_data_t **data);

void
spamm_hashed_delete_node (struct spamm_hashed_node_t **node);

unsigned int
spamm_index_row_major (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N);

unsigned int
spamm_index_column_major (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N);

unsigned int
spamm_index_column_major_2 (const unsigned int number_dimensions,
    const unsigned int N_block,
    const unsigned int *const N_lower,
    const unsigned int *const i);

unsigned int
spamm_index_norm (const unsigned int i, const unsigned int j);

unsigned int
spamm_index_kernel_block (const unsigned int i, const unsigned int j, const enum spamm_layout_t layout);

unsigned int
spamm_index_kernel_block_transpose (const unsigned int i, const unsigned int j, const enum spamm_layout_t layout);

//inline unsigned int
//spamm_dense_index (const unsigned int i, const unsigned int j)
//{
//  if (i >= SPAMM_N_KERNEL || i >= SPAMM_N_KERNEL)
//  {
//    fprintf(stderr, "i or j out of bounds\n");
//    exit(1);
//  }
//
//  return i*SPAMM_N_KERNEL+j;
//}

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

unsigned int
spamm_get_number_of_rows (const struct spamm_hashed_t *const A);

unsigned int
spamm_get_number_of_columns (const struct spamm_hashed_t *const A);

float
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

unsigned int
spamm_memory (const struct spamm_hashed_t *A);

void
spamm_hashed_multiply (const float tolerance,
    const float alpha, struct spamm_hashed_t *A, struct spamm_hashed_t *B,
    const float beta, struct spamm_hashed_t *C,
    struct spamm_timer_t *timer,
    const enum spamm_kernel_t kernel);

void
spamm_recursive_multiply_scalar (const float alpha,
    struct spamm_recursive_node_t *A,
    const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int contiguous_tier,
    const short use_linear_tree);

void
spamm_multiply (const float tolerance,
    const float alpha,
    struct spamm_matrix_t *A,
    struct spamm_matrix_t *B,
    const float beta,
    struct spamm_matrix_t *C,
    struct spamm_timer_t *timer,
    sgemm_func sgemm,
    const enum spamm_kernel_t kernel,
    unsigned int *number_products);

void
spamm_add (const float alpha,
    struct spamm_matrix_t *A,
    const float beta,
    struct spamm_matrix_t *B);

unsigned int
spamm_number_nonzero (const struct spamm_matrix_t *A);

void
spamm_print_info (const struct spamm_matrix_t *const A);

unsigned int
spamm_get_tree_depth (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int N_block,
    const short use_linear_tree);

struct spamm_matrix_t *
spamm_new (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int contiguous_tier,
    const unsigned int N_block,
    const short use_linear_tree);

struct spamm_hashed_t *
spamm_hashed_new (const unsigned int tier,
    const unsigned int kernel_tier,
    const unsigned int depth,
    const unsigned int M_lower,
    const unsigned int M_upper,
    const unsigned int N_lower,
    const unsigned int N_upper);

struct spamm_hashed_data_t *
spamm_hashed_new_data (const unsigned int tier, const unsigned int index_2D, const enum spamm_layout_t layout);

struct spamm_hashed_node_t *
spamm_hashed_new_node (const unsigned int tier, const unsigned int index_2D);

struct spamm_recursive_node_t *
spamm_recursive_new_node ();

void
spamm_print (const struct spamm_matrix_t *A);

void
spamm_print_tree (const struct spamm_matrix_t *A);

void
spamm_print_dense (const unsigned int M, const unsigned int N,
    const enum spamm_layout_t type, const float *A);

void
spamm_set (const unsigned int *const i, const float Aij, struct spamm_matrix_t *A);

void
spamm_recursive_set (const unsigned int number_dimensions,
    const unsigned int *const i,
    const unsigned int *const N,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    const unsigned int tier,
    const unsigned int contiguous_tier,
    const unsigned int kernel_tier,
    const unsigned int N_block,
    const short use_linear_tree,
    const unsigned int depth,
    const float Aij,
    struct spamm_recursive_node_t **node);

void
spamm_hashed_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_hashed_t *A);

void
spamm_uint_to_bin_string (const unsigned int width, const unsigned int i, char *result);

char *
spamm_version ();

void
spamm_construct_tree (struct spamm_hashed_t *A);

void
spamm_hashed_norm_update (struct spamm_hashed_t *A);

struct spamm_matrix_t *
spamm_convert_dense_to_spamm (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int contiguous_tier,
    const unsigned int N_block,
    const short use_linear_tree,
    const enum spamm_layout_t dense_type,
    const float *const A_dense);

struct spamm_hashed_t *
spamm_hashed_convert_dense_to_spamm (const unsigned int M, const unsigned int N,
    const enum spamm_layout_t dense_type, const float *const A_dense,
    const enum spamm_layout_t spamm_layout);

void spamm_sgemm (char * transA, char * transB,
    int *M, int *N, int *K,
    float *alpha, float *A, int *LDA, float *B, int *LDB,
    float *beta, float *C, int *LDC);

#endif
