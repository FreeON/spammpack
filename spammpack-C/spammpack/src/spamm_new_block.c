#include "spamm.h"

#include <stdlib.h>

struct spamm_data_t *
spamm_new_block (const unsigned int tier, const unsigned int index)
{
  struct spamm_data_t *data = (struct spamm_data_t*) malloc(sizeof(struct spamm_data_t));

  data->tier = tier;
  data->index = index;

  return data;
}
