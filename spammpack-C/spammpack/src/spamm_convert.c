#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//#define EXTRA_DEBUG
#define SPAMM_SET_NO_ZERO

/** Convert a dense matrix to SpAMM.
 *
 * @param M The number of rows.
 * @param N The number of columns.
 * @param dense_type The storage type of the dense matrix.
 * @param A_dense The dense matrix.
 * @param spamm_layout The layout of the SpAMM data nodes.
 *
 * @return The SpAMM matrix.
 */
struct spamm_hashed_t *
spamm_convert_dense_to_spamm (const unsigned int M, const unsigned int N,
    const enum spamm_layout_t dense_type, float *A_dense,
    const enum spamm_layout_t spamm_layout)
{
  struct spamm_hashed_t *A = NULL;
  unsigned int index;
  unsigned int i;
  unsigned int j;
  unsigned int i_blocked;
  unsigned int j_blocked;
  unsigned int i_basic;
  unsigned int j_basic;
  unsigned int i_kernel;
  unsigned int j_kernel;
  unsigned int i_tier;
  unsigned int j_tier;
  unsigned int norm_offset;
  unsigned int data_offset;
  unsigned int data_offset_transpose;
  float norm_A11;
  float norm_A12;
  float norm_A21;
  float norm_A22;
  struct spamm_hashtable_t *node_hashtable;
  struct spamm_data_t *data;
  float Aij;
  float norm2;
#ifdef SPAMM_SET_NO_ZERO
  short found_nonzero_element;
#endif

  assert(A_dense != NULL);

#ifdef EXTRA_DEBUG
  printf("creating new SpAMM %ux%u matrix\n", M, N);
#endif
  A = spamm_hashed_new(M, N, spamm_layout);

  /* Get hash table at this tier. */
  node_hashtable = A->tier_hashtable[A->kernel_tier];

  /* Store matrix elements on kernel tier. */
  for (i = 0; i < M; i += SPAMM_N_KERNEL) {
    for (j = 0; j < N; j += SPAMM_N_KERNEL)
    {
      /* Calculate the matrix block indices. */
      i_tier = i/SPAMM_N_KERNEL;
      j_tier = j/SPAMM_N_KERNEL;

      /* Construct linear index of the node on this tier. */
      index = spamm_index_2D(i_tier, j_tier);

#ifdef SPAMM_SET_NO_ZERO
      /* Check whether at least one elements is non-zero for this kernel tier
       * block. */
#ifdef EXTRA_DEBUG
      printf("[convert] analyzing kernel tier block A(%u,%u), index = %u\n", i, j, index);
#endif
      found_nonzero_element = 0;
      for (i_kernel = 0; i_kernel < SPAMM_N_KERNEL && i+i_kernel < M; i_kernel++) {
        for (j_kernel = 0; j_kernel < SPAMM_N_KERNEL && j+j_kernel < N; j_kernel++)
        {
          switch(dense_type)
          {
            case row_major:
#ifdef EXTRA_DEBUG
              printf("[convert] storing A[%u,%u] = %e\n", i+i_kernel, j+j_kernel, A_dense[spamm_index_row_major(i+i_kernel, j+j_kernel, M, N)]);
#endif
              if (A_dense[spamm_index_row_major(i+i_kernel, j+j_kernel, M, N)] != 0.0)
              {
                found_nonzero_element = 1;
              }
              break;

            case column_major:
#ifdef EXTRA_DEBUG
              printf("[convert] storing A[%u,%u] = %e\n", i+i_kernel, j+j_kernel, A_dense[spamm_index_column_major(i+i_kernel, j+j_kernel, M, N)]);
#endif
              if (A_dense[spamm_index_column_major(i+i_kernel, j+j_kernel, M, N)] != 0.0)
              {
                found_nonzero_element = 1;
              }
              break;

            default:
              printf("unknown type\n");
              exit(1);
              break;
          }

          if (found_nonzero_element == 1) { break; }
        }
        if (found_nonzero_element == 1) { break; }
      }

      if (found_nonzero_element == 0)
      {
#ifdef EXTRA_DEBUG
        printf("[convert] all elements are zero in this block\n");
#endif
        continue;
      }

#ifdef EXTRA_DEBUG
      else
      {
        printf("[convert] found at least 1 non zero element\n");
      }
#endif

#endif

      if ((data = spamm_hashtable_lookup(node_hashtable, index)) == NULL)
      {
        data = spamm_new_block(A->kernel_tier, index, A->layout);
        spamm_hashtable_insert(node_hashtable, index, data);
      }

      /* Set all elements in this kernel block. */
      for (i_blocked = 0; i_blocked < SPAMM_N_KERNEL_BLOCKED && i+SPAMM_N_BLOCK*i_blocked < M; i_blocked++) {
        for (j_blocked = 0; j_blocked < SPAMM_N_KERNEL_BLOCKED && j+SPAMM_N_BLOCK*j_blocked < N; j_blocked++)
        {
          norm_offset = spamm_index_norm(i_blocked, j_blocked);

          /* Reset norm2. */
          norm2 = 0.0;

          for (i_basic = 0; i_basic < SPAMM_N_BLOCK && i+SPAMM_N_BLOCK*i_blocked+i_basic < M; i_basic++) {
            for (j_basic = 0; j_basic < SPAMM_N_BLOCK && j+SPAMM_N_BLOCK*j_blocked+j_basic < N; j_basic++)
            {
              /* Calculate offsets into the matrix data. */
              data_offset = spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, A->layout);
              data_offset_transpose = spamm_index_kernel_block_transpose_hierarchical(i_blocked, j_blocked, i_basic, j_basic, A->layout);

              /* Get matrix elements from dense matrix. */
              switch(dense_type)
              {
                case row_major:
                  Aij = A_dense[spamm_index_row_major(i+SPAMM_N_BLOCK*i_blocked+i_basic, j+SPAMM_N_BLOCK*j_blocked+j_basic, M, N)];
                  break;

                case column_major:
                  Aij = A_dense[spamm_index_column_major(i+SPAMM_N_BLOCK*i_blocked+i_basic, j+SPAMM_N_BLOCK*j_blocked+j_basic, M, N)];
                  break;

                default:
                  printf("unknown type\n");
                  exit(1);
              }

              /* Set new value. */
              data->block_dense[data_offset] = Aij;
              data->block_dense_store[data_offset] = Aij;
              data->block_dense_transpose[data_offset_transpose] = Aij;

              data->block_dense_dilated[4*data_offset+0] = Aij;
              data->block_dense_dilated[4*data_offset+1] = Aij;
              data->block_dense_dilated[4*data_offset+2] = Aij;
              data->block_dense_dilated[4*data_offset+3] = Aij;

              /* Update norm. */
              norm2 += Aij*Aij;
            }
          }
          data->norm2[norm_offset] = norm2;
          data->norm[norm_offset] = sqrt(norm2);
        }
      }

      /* Loop over upper tier. */
      norm_A11 = sqrt(
          data->norm2[spamm_index_norm(0, 0)]+
          data->norm2[spamm_index_norm(0, 1)]+
          data->norm2[spamm_index_norm(1, 0)]+
          data->norm2[spamm_index_norm(1, 1)]);
      norm_A12 = sqrt(
          data->norm2[spamm_index_norm(0, 2)]+
          data->norm2[spamm_index_norm(0, 3)]+
          data->norm2[spamm_index_norm(1, 2)]+
          data->norm2[spamm_index_norm(1, 3)]);
      norm_A21 = sqrt(
          data->norm2[spamm_index_norm(2, 0)]+
          data->norm2[spamm_index_norm(2, 1)]+
          data->norm2[spamm_index_norm(3, 0)]+
          data->norm2[spamm_index_norm(3, 1)]);
      norm_A22 = sqrt(
          data->norm2[spamm_index_norm(2, 2)]+
          data->norm2[spamm_index_norm(2, 3)]+
          data->norm2[spamm_index_norm(3, 2)]+
          data->norm2[spamm_index_norm(3, 3)]);

      data->norm_upper[0] = norm_A11;
      data->norm_upper[1] = norm_A12;
      data->norm_upper[2] = norm_A11;
      data->norm_upper[3] = norm_A12;
      data->norm_upper[4] = norm_A21;
      data->norm_upper[5] = norm_A22;
      data->norm_upper[6] = norm_A21;
      data->norm_upper[7] = norm_A22;

      data->norm_upper_transpose[0] = norm_A11;
      data->norm_upper_transpose[1] = norm_A21;
      data->norm_upper_transpose[2] = norm_A12;
      data->norm_upper_transpose[3] = norm_A22;
      data->norm_upper_transpose[4] = norm_A11;
      data->norm_upper_transpose[5] = norm_A21;
      data->norm_upper_transpose[6] = norm_A12;
      data->norm_upper_transpose[7] = norm_A22;

      /* Update node norm. */
      for (i_blocked = 0; i_blocked < SPAMM_N_KERNEL_BLOCKED; i_blocked++) {
        for (j_blocked = 0; j_blocked < SPAMM_N_KERNEL_BLOCKED; j_blocked++)
        {
          data->node_norm2 += data->norm2[spamm_index_norm(i_blocked, j_blocked)];
        }
      }
      data->node_norm = sqrt(data->node_norm2);
    }
  }

  /* Construct the rest of the tree. */
  spamm_construct_tree(A);

  return A;
}
