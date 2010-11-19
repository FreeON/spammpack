#include "spamm.h"

#include <stdlib.h>

void
spamm_delete_node (struct spamm_node_t **node)
{
  free(*node);
  *node = NULL;
}
