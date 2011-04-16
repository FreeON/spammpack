/** @file */

#ifndef __SPAMM_H
#define __SPAMM_H

/* Include some configuration options. */
#include "spamm_config.h"
#include "spamm_error.h"
#include "spamm_hashtable.h"
#include "spamm_introspection.h"
#include "spamm_list.h"
#include "spamm_kernel.h"
#include "spamm_timer.h"
#include "spamm_convert.h"
#include "spamm_types.h"

void *
spamm_allocate (size_t size);

int
spamm_check (const struct spamm_t *A);

void
spamm_delete (struct spamm_t **A);

void
spamm_delete_block (struct spamm_data_t **data);

void
spamm_delete_node (struct spamm_node_t **node);

unsigned int
spamm_index_kernel_block (const unsigned int i, const unsigned int j, const enum spamm_layout_t layout);

unsigned int
spamm_index_row_major (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N);

unsigned int
spamm_index_column_major (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N);

unsigned int
spamm_index_norm (const unsigned int i, const unsigned int j);

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
spamm_index_kernel_block_hierarchical_1 (const unsigned int i_block,
    const unsigned int j_block, const unsigned int i,
    const unsigned int j, const enum spamm_layout_t layout);

float
spamm_get (const unsigned int i, const unsigned int j, const struct spamm_t *A);

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
spamm_memory (const struct spamm_t *A);

void
spamm_multiply (const float tolerance,
    const float alpha, struct spamm_t *A, struct spamm_t *B,
    const float beta, struct spamm_t *C,
    const enum spamm_timer_type_t timer_type,
    const enum spamm_kernel_t kernel);

unsigned int
spamm_number_nonzero (const struct spamm_t *A);

struct spamm_t *
spamm_new (const unsigned int M, const unsigned int N, const enum spamm_layout_t layout);

struct spamm_data_t *
spamm_new_block (const unsigned int tier, const unsigned int index_2D, const enum spamm_layout_t layout);

struct spamm_node_t *
spamm_new_node (const unsigned int tier, const unsigned int index_2D);

void
spamm_print (const struct spamm_t *A);

void
spamm_print_tree (const struct spamm_t *A);

void
spamm_print_dense (const unsigned int M, const unsigned int N,
    const enum spamm_layout_t type, const float *A);

void
spamm_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_t *A);

void
spamm_uint_to_bin_string (const unsigned int width, const unsigned int i, char *result);

char *
spamm_version ();

void
spamm_prune (struct spamm_t *A);

void
spamm_expand (struct spamm_t *A);

void
spamm_construct_tree (struct spamm_t *A);

#endif
