#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Convert a dense matrix to SpAMM.
 *
 * @param number_dimensions The number of dimensions.
 * @param N The number of rows/columns.
 * @param N_linear The linear tier (see spamm_new()).
 * @param contiguous_tier The contiguous tier (see spamm_new()).
 * @param dense_type The storage type of the dense matrix.
 * @param A_dense The dense matrix.
 * @param spamm_layout The layout of the SpAMM data nodes.
 *
 * @return The SpAMM matrix.
 */
struct spamm_matrix_t *
spamm_convert_dense_to_spamm (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int contiguous_tier,
    const short use_linear_tree,
    const enum spamm_layout_t dense_type,
    const float *const A_dense,
    const enum spamm_layout_t spamm_layout)
{
  struct spamm_matrix_t *A = NULL;
  unsigned int *i;

  assert(A_dense != NULL);

  if (number_dimensions != 2)
  {
    SPAMM_FATAL("can not handle this case\n");
  }

  A = spamm_new(number_dimensions, N, contiguous_tier, use_linear_tree, spamm_layout);

  i = calloc(number_dimensions, sizeof(unsigned int));
  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
    {
      switch(dense_type)
      {
        case row_major:
          spamm_set(i, A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])], A);
          break;

        case column_major:
          spamm_set(i, A_dense[spamm_index_column_major(i[0], i[1], N[0], N[1])], A);
          break;

        default:
          SPAMM_FATAL("unknown type\n");
          break;
      }
    }
  }
  free(i);

  return A;
}

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
//struct spamm_hashed_t *
//spamm_hashed_convert_dense_to_spamm (const unsigned int M, const unsigned int N,
//    const enum spamm_layout_t dense_type, const float *const A_dense,
//    const enum spamm_layout_t spamm_layout)
//{
//  struct spamm_hashed_t *A = NULL;
//  unsigned int index;
//  unsigned int i;
//  unsigned int j;
//  unsigned int i_blocked;
//  unsigned int j_blocked;
//  unsigned int i_basic;
//  unsigned int j_basic;
//  unsigned int i_kernel;
//  unsigned int j_kernel;
//  unsigned int i_tier;
//  unsigned int j_tier;
//  unsigned int norm_offset;
//  unsigned int data_offset;
//  unsigned int data_offset_transpose;
//  struct spamm_hashtable_t *node_hashtable;
//  struct spamm_hashed_data_t *data;
//  float Aij;
//  float norm2;
//#ifdef SPAMM_SET_NO_ZERO
//  short found_nonzero_element;
//#endif
//
//  assert(A_dense != NULL);
//
//#ifdef EXTRA_DEBUG
//  printf("creating new SpAMM %ux%u matrix\n", M, N);
//#endif
//  A = spamm_hashed_new(M, N, spamm_layout);
//
//  /* Get hash table at this tier. */
//  node_hashtable = A->tier_hashtable[A->kernel_tier];
//
//  /* Store matrix elements on kernel tier. */
//  for (i = 0; i < M; i += SPAMM_N_KERNEL) {
//    for (j = 0; j < N; j += SPAMM_N_KERNEL)
//    {
//      /* Calculate the matrix block indices. */
//      i_tier = i/SPAMM_N_KERNEL;
//      j_tier = j/SPAMM_N_KERNEL;
//
//      /* Construct linear index of the node on this tier. */
//      index = spamm_index_2D(i_tier, j_tier);
//
//#ifdef SPAMM_SET_NO_ZERO
//      /* Check whether at least one elements is non-zero for this kernel tier
//       * block. */
//#ifdef EXTRA_DEBUG
//      printf("[convert] analyzing kernel tier block A(%u,%u), index = %u\n", i, j, index);
//#endif
//      found_nonzero_element = 0;
//      for (i_kernel = 0; i_kernel < SPAMM_N_KERNEL && i+i_kernel < M; i_kernel++) {
//        for (j_kernel = 0; j_kernel < SPAMM_N_KERNEL && j+j_kernel < N; j_kernel++)
//        {
//          switch(dense_type)
//          {
//            case row_major:
//#ifdef EXTRA_DEBUG
//              printf("[convert] storing A[%u,%u] = %e\n", i+i_kernel, j+j_kernel, A_dense[spamm_index_row_major(i+i_kernel, j+j_kernel, M, N)]);
//#endif
//              if (A_dense[spamm_index_row_major(i+i_kernel, j+j_kernel, M, N)] != 0.0)
//              {
//                found_nonzero_element = 1;
//              }
//              break;
//
//            case column_major:
//#ifdef EXTRA_DEBUG
//              printf("[convert] storing A[%u,%u] = %e\n", i+i_kernel, j+j_kernel, A_dense[spamm_index_column_major(i+i_kernel, j+j_kernel, M, N)]);
//#endif
//              if (A_dense[spamm_index_column_major(i+i_kernel, j+j_kernel, M, N)] != 0.0)
//              {
//                found_nonzero_element = 1;
//              }
//              break;
//
//            default:
//              SPAMM_FATAL("unknown type\n");
//              break;
//          }
//
//          if (found_nonzero_element == 1) { break; }
//        }
//        if (found_nonzero_element == 1) { break; }
//      }
//
//      if (found_nonzero_element == 0)
//      {
//#ifdef EXTRA_DEBUG
//        printf("[convert] all elements are zero in this block\n");
//#endif
//        continue;
//      }
//
//#ifdef EXTRA_DEBUG
//      else
//      {
//        printf("[convert] found at least 1 non zero element\n");
//      }
//#endif
//
//#endif
//
//      if ((data = spamm_hashtable_lookup(node_hashtable, index)) == NULL)
//      {
//        data = spamm_hashed_new_data(A->kernel_tier, index, A->layout);
//        spamm_hashtable_insert(node_hashtable, index, data);
//      }
//
//      /* Set all elements in this kernel block. */
//      for (i_blocked = 0; i_blocked < SPAMM_N_KERNEL_BLOCKED && i+SPAMM_N_BLOCK*i_blocked < M; i_blocked++) {
//        for (j_blocked = 0; j_blocked < SPAMM_N_KERNEL_BLOCKED && j+SPAMM_N_BLOCK*j_blocked < N; j_blocked++)
//        {
//          norm_offset = spamm_index_norm(i_blocked, j_blocked);
//
//          /* Reset norm2. */
//          norm2 = 0.0;
//
//          for (i_basic = 0; i_basic < SPAMM_N_BLOCK && i+SPAMM_N_BLOCK*i_blocked+i_basic < M; i_basic++) {
//            for (j_basic = 0; j_basic < SPAMM_N_BLOCK && j+SPAMM_N_BLOCK*j_blocked+j_basic < N; j_basic++)
//            {
//              /* Calculate offsets into the matrix data. */
//              data_offset = spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, A->layout);
//              data_offset_transpose = spamm_index_kernel_block_transpose_hierarchical(i_blocked, j_blocked, i_basic, j_basic, A->layout);
//
//              /* Get matrix elements from dense matrix. */
//              switch(dense_type)
//              {
//                case row_major:
//                  Aij = A_dense[spamm_index_row_major(i+SPAMM_N_BLOCK*i_blocked+i_basic, j+SPAMM_N_BLOCK*j_blocked+j_basic, M, N)];
//                  break;
//
//                case column_major:
//                  Aij = A_dense[spamm_index_column_major(i+SPAMM_N_BLOCK*i_blocked+i_basic, j+SPAMM_N_BLOCK*j_blocked+j_basic, M, N)];
//                  break;
//
//                default:
//                  SPAMM_FATAL("unknown type\n");
//                  break;
//              }
//
//              /* Set new value. */
//              data->block_dense[data_offset] = Aij;
//              data->block_dense_store[data_offset] = Aij;
//              data->block_dense_transpose[data_offset_transpose] = Aij;
//
//              data->block_dense_dilated[4*data_offset+0] = Aij;
//              data->block_dense_dilated[4*data_offset+1] = Aij;
//              data->block_dense_dilated[4*data_offset+2] = Aij;
//              data->block_dense_dilated[4*data_offset+3] = Aij;
//
//              /* Update norm. */
//              norm2 += Aij*Aij;
//            }
//          }
//          data->norm2[norm_offset] = norm2;
//          data->norm[norm_offset] = sqrt(norm2);
//        }
//      }
//
//      /* Update node norm. */
//      for (i_blocked = 0; i_blocked < SPAMM_N_KERNEL_BLOCKED; i_blocked++) {
//        for (j_blocked = 0; j_blocked < SPAMM_N_KERNEL_BLOCKED; j_blocked++)
//        {
//          data->node_norm2 += data->norm2[spamm_index_norm(i_blocked, j_blocked)];
//        }
//      }
//      data->node_norm = sqrt(data->node_norm2);
//    }
//  }
//
//  /* Construct the rest of the tree. */
//  spamm_construct_tree(A);
//
//  return A;
//}
