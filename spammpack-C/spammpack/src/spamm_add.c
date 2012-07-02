/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Add two spamm matrices. @f$ A \leftarrow \alpha A + \beta B @f$.
 *
 * @param alpha The factor @f$ \alpha @f$.
 * @param A The matrix A.
 * @param beta The factor @f$ \beta @f$.
 * @param B The matrix B.
 */
void
spamm_hashed_add (const float alpha,
    struct spamm_hashed_t *A,
    const float beta,
    struct spamm_hashed_t *B)
{
  struct spamm_list_t *A_block_keys;
  struct spamm_list_t *B_block_keys;

  struct spamm_node_t *node;
  struct spamm_hashed_data_t *A_data;
  struct spamm_hashed_data_t *B_data;

  unsigned int A_index;
  unsigned int B_index;

  int i, j, k;
  int i_block, j_block;

  assert(A != NULL);
  assert(B != NULL);

  if (A->layout != B->layout)
  {
    printf("[add] inconsisten layout in matrices\n");
    exit(1);
  }

  if (A->M != B->M)
  {
    printf("[add] mismatch of number of rows\n");
    exit(1);
  }

  if (A->N != B->N)
  {
    printf("[add] mismatch of number of columns\n");
    exit(1);
  }

  /* Print out some information. */
  printf("[add] alpha = %e, beta = %e\n", alpha, beta);

  /* Go to lowest tier and start adding submatrix blocks. */
  A_block_keys = spamm_hashtable_keys(A->tier_hashtable[A->kernel_tier]);
  B_block_keys = spamm_hashtable_keys(B->tier_hashtable[A->kernel_tier]);

  /* Sort the keys. */
  spamm_list_sort_index(A_block_keys, spamm_list_compare_int);
  spamm_list_sort_index(B_block_keys, spamm_list_compare_int);

  for (i = 0, j = 0; ; )
  {
    if (i < spamm_list_length(A_block_keys))
    {
      A_index = spamm_list_get_index(A_block_keys, i);
    }

    else
    {
      A_index = UINT_MAX;
    }

    if (j < spamm_list_length(B_block_keys))
    {
      B_index = spamm_list_get_index(B_block_keys, j);
    }

    else
    {
      B_index = UINT_MAX;
    }

    if (i >= spamm_list_length(A_block_keys) && j >= spamm_list_length(B_block_keys))
    {
      /* Done. */
      break;
    }

    if (A_index < B_index)
    {
      A_data = spamm_hashtable_lookup(A->tier_hashtable[A->kernel_tier], A_index);
      for (k = 0; k < SPAMM_N_KERNEL*SPAMM_N_KERNEL; k++)
      {
        A_data->block_dense[k] *= alpha;
      }
      A_data->node_norm2 *= alpha*alpha;
      A_data->node_norm *= alpha;

      i++;

      /* Update block norms. */
      for (i_block = 0; i_block < SPAMM_N_KERNEL_BLOCKED; i_block++) {
        for (j_block = 0; j_block < SPAMM_N_KERNEL_BLOCKED; j_block++)
        {
          A_data->norm2[spamm_index_norm(i_block, j_block)] = 0;
          for (k = 0; k < SPAMM_N_BLOCK*SPAMM_N_BLOCK; k++)
          {
            A_data->norm2[spamm_index_norm(i_block, j_block)] +=
              A_data->block_dense[spamm_index_kernel_block(i_block, j_block, A->layout)+k]
              *A_data->block_dense[spamm_index_kernel_block(i_block, j_block, A->layout)+k];
          }
          A_data->norm[spamm_index_norm(i_block, j_block)] = sqrt(A_data->norm[spamm_index_norm(i_block, j_block)]);
        }
      }
    }

    else if (B_index < A_index)
    {
      B_data = spamm_hashtable_lookup(B->tier_hashtable[B->kernel_tier], spamm_list_get_index(B_block_keys, j));

      /* Create new node. */
      A_data = spamm_hashed_new_data(A->kernel_tier, B_index, A->layout);
      A_data->node_norm2 = 0;
      for (k = 0; k < SPAMM_N_KERNEL*SPAMM_N_KERNEL; k++)
      {
        A_data->block_dense[k] = beta*B_data->block_dense[k];
        A_data->node_norm2 += A_data->block_dense[k]*A_data->block_dense[k];
      }
      A_data->node_norm = sqrt(A_data->node_norm2);

      j++;

      /* Update block norms. */
      for (i_block = 0; i_block < SPAMM_N_KERNEL_BLOCKED; i_block++) {
        for (j_block = 0; j_block < SPAMM_N_KERNEL_BLOCKED; j_block++)
        {
          A_data->norm2[spamm_index_norm(i_block, j_block)] = 0;
          for (k = 0; k < SPAMM_N_BLOCK*SPAMM_N_BLOCK; k++)
          {
            A_data->norm2[spamm_index_norm(i_block, j_block)] +=
              A_data->block_dense[spamm_index_kernel_block(i_block, j_block, A->layout)+k]
              *A_data->block_dense[spamm_index_kernel_block(i_block, j_block, A->layout)+k];
          }
          A_data->norm[spamm_index_norm(i_block, j_block)] = sqrt(A_data->norm[spamm_index_norm(i_block, j_block)]);
        }
      }

      /* Create new block in A and store it. */
      spamm_hashtable_insert(A->tier_hashtable[A->kernel_tier], B_index, A_data);
    }

    else if (A_index == B_index)
    {
      A_data = spamm_hashtable_lookup(A->tier_hashtable[A->kernel_tier], spamm_list_get_index(A_block_keys, i));
      B_data = spamm_hashtable_lookup(B->tier_hashtable[B->kernel_tier], spamm_list_get_index(B_block_keys, j));

      /* Add block data. */
      A_data->node_norm2 = 0;
      for (k = 0; k < SPAMM_N_KERNEL*SPAMM_N_KERNEL; k++)
      {
        A_data->block_dense[k] = alpha*A_data->block_dense[k]+beta*B_data->block_dense[k];
        A_data->node_norm2 += A_data->block_dense[k]*A_data->block_dense[k];
      }
      A_data->node_norm = sqrt(A_data->node_norm2);

      i++;
      j++;

      /* Update block norms. */
      for (i_block = 0; i_block < SPAMM_N_KERNEL_BLOCKED; i_block++) {
        for (j_block = 0; j_block < SPAMM_N_KERNEL_BLOCKED; j_block++)
        {
          A_data->norm2[spamm_index_norm(i_block, j_block)] = 0;
          for (k = 0; k < SPAMM_N_BLOCK*SPAMM_N_BLOCK; k++)
          {
            A_data->norm2[spamm_index_norm(i_block, j_block)] +=
              A_data->block_dense[spamm_index_kernel_block(i_block, j_block, A->layout)+k]
              *A_data->block_dense[spamm_index_kernel_block(i_block, j_block, A->layout)+k];
          }
          A_data->norm[spamm_index_norm(i_block, j_block)] = sqrt(A_data->norm[spamm_index_norm(i_block, j_block)]);
        }
      }
    }

    else
    {
      printf("[%s:%i] I should not be here\n", __FILE__, __LINE__);
      exit(1);
    }
  }

  /* Insert into upper tier hashtables. */
}
