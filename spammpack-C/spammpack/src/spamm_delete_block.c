#include "spamm.h"
#include <stdlib.h>

void
spamm_delete_block (struct spamm_data_t **data)
{
  free(*data);
  *data = NULL;
}
