#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

void
spamm_delete_node (struct spamm_node_t *node)
{
  int i;

  for (i = 0; i < node->child_M*node->child_N; ++i)
  {
  }
}

void
spamm_delete (struct spamm_t *A)
{
  assert(A != NULL);

  /* Recurse down and free() nodes. */
  spamm_delete_node(A->root);
}
