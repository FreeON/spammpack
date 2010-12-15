#include "spamm.h"
#include <stdlib.h>

/** Allocate a new data node of a matrix tree.
 *
 * @param tier The tier this node will be on.
 * @param index_2D The 2D linear matrix index of this node.
 * @param index_3D_ik0 The 3D convolution index in ik0 order.
 * @param index_3D_0kj The 3D convolution index in 0kj order.
 *
 * @return A pointer to the newly allocated node.
 */
struct spamm_data_t *
spamm_new_block (const unsigned int tier, const unsigned int index_2D,
    const unsigned int index_3D_ik0, const unsigned int index_3D_0kj)
{
  int i, j;
  struct spamm_data_t *data = (struct spamm_data_t*) malloc(sizeof(struct spamm_data_t));

  data->tier = tier;
  data->index_2D = index_2D;
  data->index_3D_ik0 = index_3D_ik0;
  data->index_3D_0kj = index_3D_0kj;

  /* Set matrix elements to zero. */
  for (i = 0; i < SPAMM_N_KERNEL; i++) {
    for (j = 0; j < SPAMM_N_KERNEL; j++)
    {
      data->block_dense[i*SPAMM_N_KERNEL+j] = 0.0;

      data->block_dense_dilated[4*(i*SPAMM_N_KERNEL+j)+0] = 0.0;
      data->block_dense_dilated[4*(i*SPAMM_N_KERNEL+j)+1] = 0.0;
      data->block_dense_dilated[4*(i*SPAMM_N_KERNEL+j)+2] = 0.0;
      data->block_dense_dilated[4*(i*SPAMM_N_KERNEL+j)+3] = 0.0;
    }
  }

  for (i = 0; i < SPAMM_N_KERNEL_BLOCK*SPAMM_N_KERNEL_BLOCK; i++)
  {
    data->norm[i] = 0.0;
    data->norm2[i] = 0.0;
  }

  data->node_norm = 0.0;
  data->node_norm2 = 0.0;

  return data;
}
