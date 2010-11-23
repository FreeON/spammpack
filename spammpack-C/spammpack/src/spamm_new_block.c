#include "spamm.h"
#include <stdlib.h>

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

  for (i = 0; i < SPAMM_N_KERNEL; i++) {
    for (j = 0; j < SPAMM_N_KERNEL; j++)
    {
      data->block_dense[i*SPAMM_N_KERNEL+j] = 0.0;
    }
  }

  for (i = 0; i < SPAMM_N_KERNEL_BLOCK*SPAMM_N_KERNEL_BLOCK; i++)
  {
    data->norm[i] = 0.0;
  }

  return data;
}
