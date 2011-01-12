#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Convert a dense matrix to SpAMM.
 *
 * @param M The number of rows.
 * @param N The number of columns.
 * @param A_dense The dense matrix.
 *
 * @return The SpAMM matrix.
 */
  struct spamm_t *
spamm_convert_dense_to_spamm (const unsigned int M, const unsigned int N,
    const enum spamm_dense_type_t type, float *A_dense)
{
  struct spamm_t *A = NULL;
  unsigned int tier;
  unsigned int next_tier;
  unsigned int reverse_tier;
  unsigned int index;
  unsigned int parent_index;
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
  struct spamm_hashtable_t *node_hashtable;
  struct spamm_hashtable_t *next_tier_hashtable;
  struct spamm_node_t *node;
  struct spamm_node_t *parent_node;
  struct spamm_data_t *data;
  struct spamm_list_t *tier_indices;
  float Aij;
  float norm2;

  assert(A_dense != NULL);

  A = spamm_new(M, N);

  /* Get hash table at this tier. */
  node_hashtable = A->tier_hashtable[A->kernel_tier];

  for (i = 0; i < M; i += SPAMM_N_KERNEL) {
    for (j = 0; j < N; j += SPAMM_N_KERNEL)
    {
      /* Calculate the matrix block indices. */
      i_tier = i/SPAMM_N_KERNEL;
      j_tier = j/SPAMM_N_KERNEL;

      /* Construct linear index of the node on this tier. */
      index = spamm_index_2D(i_tier, j_tier);

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
          norm2 = 0.0;

          for (i_block = 0; i_block < SPAMM_N_BLOCK && i+SPAMM_N_BLOCK*i_kernel+i_block < M; i_block++) {
            for (j_block = 0; j_block < SPAMM_N_BLOCK && j+SPAMM_N_BLOCK*j_kernel+j_block < N; j_block++)
            {
              /* Calculate offsets into the matrix data. */
              data_offset = spamm_index_kernel_block(SPAMM_N_BLOCK*i_kernel+i_block, SPAMM_N_BLOCK*j_kernel+j_block);

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

  for (tier = 0; tier <= A->kernel_tier-1; tier++)
  {
    reverse_tier = A->kernel_tier-tier;
    node_hashtable = A->tier_hashtable[reverse_tier];
    next_tier = reverse_tier-1;
    next_tier_hashtable = A->tier_hashtable[next_tier];

    tier_indices = spamm_hashtable_keys(node_hashtable);

    if (reverse_tier == A->kernel_tier)
    {
      for (i = 0; i < spamm_list_length(tier_indices); i++)
      {
        /* Get block. */
        data = spamm_hashtable_lookup(node_hashtable, spamm_list_get(tier_indices, i));

        /* Construct index of parent node. */
        parent_index = data->index_2D >> 2;

        /* Get parent node. */
        parent_node = spamm_hashtable_lookup(next_tier_hashtable, parent_index);

        if (parent_node == NULL)
        {
          parent_node = spamm_new_node(next_tier, parent_index);
          spamm_hashtable_insert(next_tier_hashtable, parent_index, parent_node);
        }

        parent_node->norm2 += data->node_norm2;
      }
    }

    else
    {
      for (i = 0; i < spamm_list_length(tier_indices); i++)
      {
        /* Get block. */
        node = spamm_hashtable_lookup(node_hashtable, spamm_list_get(tier_indices, i));

        /* Construct index of parent node. */
        parent_index = node->index_2D >> 2;

        /* Get parent node. */
        parent_node = spamm_hashtable_lookup(next_tier_hashtable, parent_index);

        if (parent_node == NULL)
        {
          parent_node = spamm_new_node(next_tier, parent_index);
          spamm_hashtable_insert(next_tier_hashtable, parent_index, parent_node);
        }

        parent_node->norm2 += node->norm2;
      }
    }
    spamm_list_delete(&tier_indices);

    tier_indices = spamm_hashtable_keys(next_tier_hashtable);
    for (i = 0; i < spamm_list_length(tier_indices); i++)
    {
      /* Get block. */
      node = spamm_hashtable_lookup(next_tier_hashtable, spamm_list_get(tier_indices, i));
      node->norm = sqrt(node->norm2);
    }
    spamm_list_delete(&tier_indices);
  }

  return A;
}
