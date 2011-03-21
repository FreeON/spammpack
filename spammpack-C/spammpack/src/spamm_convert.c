#include "config.h"
#include "spamm.h"
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
 * @param type The storage type of the dense matrix.
 * @param A_dense The dense matrix.
 *
 * @return The SpAMM matrix.
 */
struct spamm_t *
spamm_convert_dense_to_spamm (const unsigned int M, const unsigned int N,
    const enum spamm_dense_type_t type, float *A_dense)
{
  struct spamm_t *A = NULL;
  unsigned int index;
  unsigned int i;
  unsigned int j;
  unsigned int i_block;
  unsigned int j_block;
  unsigned int i_kernel;
  unsigned int j_kernel;
  unsigned int i_tier;
  unsigned int j_tier;
  unsigned int norm_offset;
  unsigned int data_offset;
#ifdef SPAMM_USE_TRANSPOSE
  unsigned int data_offset_transpose;
#endif
#ifdef SPAMM_USE_HIERARCHICAL_NORM
  float norm_A11;
  float norm_A12;
  float norm_A21;
  float norm_A22;
#endif
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
  A = spamm_new(M, N);

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
          switch(type)
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
        data = spamm_new_block(A->kernel_tier, index);
        spamm_hashtable_insert(node_hashtable, index, data);
      }

      /* Set all elements in this kernel block. */
      for (i_kernel = 0; i_kernel < SPAMM_N_KERNEL_BLOCK && i+SPAMM_N_BLOCK*i_kernel < M; i_kernel++) {
        for (j_kernel = 0; j_kernel < SPAMM_N_KERNEL_BLOCK && j+SPAMM_N_BLOCK*j_kernel < N; j_kernel++)
        {
          norm_offset = spamm_index_norm(SPAMM_N_BLOCK*i_kernel, SPAMM_N_BLOCK*j_kernel);

          /* Reset norm2. */
          norm2 = 0.0;

          for (i_block = 0; i_block < SPAMM_N_BLOCK && i+SPAMM_N_BLOCK*i_kernel+i_block < M; i_block++) {
            for (j_block = 0; j_block < SPAMM_N_BLOCK && j+SPAMM_N_BLOCK*j_kernel+j_block < N; j_block++)
            {
              /* Calculate offsets into the matrix data. */
              data_offset = spamm_index_kernel_block(SPAMM_N_BLOCK*i_kernel+i_block, SPAMM_N_BLOCK*j_kernel+j_block);
#ifdef SPAMM_USE_TRANSPOSE
              data_offset_transpose = spamm_index_kernel_block(SPAMM_N_BLOCK*i_kernel+j_block, SPAMM_N_BLOCK*j_kernel+i_block);
#endif

              /* Get matrix elements from dense matrix. */
              switch(type)
              {
                case row_major:
                  Aij = A_dense[spamm_index_row_major(i+SPAMM_N_BLOCK*i_kernel+i_block, j+SPAMM_N_BLOCK*j_kernel+j_block, M, N)];
                  break;

                case column_major:
                  Aij = A_dense[spamm_index_column_major(i+SPAMM_N_BLOCK*i_kernel+i_block, j+SPAMM_N_BLOCK*j_kernel+j_block, M, N)];
                  break;

                default:
                  printf("unknown type\n");
                  exit(1);
              }

              /* Set new value. */
              data->block_dense[data_offset] = Aij;
#ifdef SPAMM_USE_TRANSPOSE
              data->block_dense_transpose[data_offset_transpose] = Aij;
#endif

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

#ifdef SPAMM_USE_HIERARCHICAL_NORM
      /* Loop over upper tier. */
      norm_A11 = sqrt(
          data->norm2[spamm_index_norm(0*SPAMM_N_BLOCK, 0*SPAMM_N_BLOCK)]+
          data->norm2[spamm_index_norm(0*SPAMM_N_BLOCK, 1*SPAMM_N_BLOCK)]+
          data->norm2[spamm_index_norm(1*SPAMM_N_BLOCK, 0*SPAMM_N_BLOCK)]+
          data->norm2[spamm_index_norm(1*SPAMM_N_BLOCK, 1*SPAMM_N_BLOCK)]);
      norm_A12 = sqrt(
          data->norm2[spamm_index_norm(0*SPAMM_N_BLOCK, 2*SPAMM_N_BLOCK)]+
          data->norm2[spamm_index_norm(0*SPAMM_N_BLOCK, 3*SPAMM_N_BLOCK)]+
          data->norm2[spamm_index_norm(1*SPAMM_N_BLOCK, 2*SPAMM_N_BLOCK)]+
          data->norm2[spamm_index_norm(1*SPAMM_N_BLOCK, 3*SPAMM_N_BLOCK)]);
      norm_A21 = sqrt(
          data->norm2[spamm_index_norm(2*SPAMM_N_BLOCK, 0*SPAMM_N_BLOCK)]+
          data->norm2[spamm_index_norm(2*SPAMM_N_BLOCK, 1*SPAMM_N_BLOCK)]+
          data->norm2[spamm_index_norm(3*SPAMM_N_BLOCK, 0*SPAMM_N_BLOCK)]+
          data->norm2[spamm_index_norm(3*SPAMM_N_BLOCK, 1*SPAMM_N_BLOCK)]);
      norm_A22 = sqrt(
          data->norm2[spamm_index_norm(2*SPAMM_N_BLOCK, 2*SPAMM_N_BLOCK)]+
          data->norm2[spamm_index_norm(2*SPAMM_N_BLOCK, 3*SPAMM_N_BLOCK)]+
          data->norm2[spamm_index_norm(3*SPAMM_N_BLOCK, 2*SPAMM_N_BLOCK)]+
          data->norm2[spamm_index_norm(3*SPAMM_N_BLOCK, 3*SPAMM_N_BLOCK)]);

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
#endif

      /* Update node norm. */
      for (i_block = 0; i_block < SPAMM_N_KERNEL_BLOCK; i_block++) {
        for (j_block = 0; j_block < SPAMM_N_KERNEL_BLOCK; j_block++)
        {
          data->node_norm2 += data->norm2[spamm_index_row_major(i_block, j_block, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)];
        }
      }
      data->node_norm = sqrt(data->node_norm2);
    }
  }

  /* Construct the rest of the tree. */
  spamm_construct_tree(A);

  return A;
}
