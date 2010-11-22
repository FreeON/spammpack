#include "spamm.h"
#include <stdlib.h>

struct spamm_data_t *
spamm_new_block (const unsigned int tier, const unsigned int index)
{
  int i, j;
  struct spamm_data_t *data = (struct spamm_data_t*) malloc(sizeof(struct spamm_data_t));

  data->tier = tier;
  data->index = index;

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
