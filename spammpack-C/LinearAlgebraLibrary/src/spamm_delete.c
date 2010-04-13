#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

void
spamm_delete_node (struct spamm_node_t *node)
{
  int i;

  for (i = 0; i < node->child_M*node->child_N; ++i)
  {
    if(node->child[i] != NULL)
    {
      spamm_delete_node(node->child[i]);
      free(node->child[i]);
    }
  }
}

void
spamm_delete (struct spamm_t *A)
{
  assert(A != NULL);

  /* Recurse down and free() nodes. */
  if (A->root != NULL)
  {
    spamm_delete_node(A->root);
  }

  free(A->root);

  A->M = 0;
  A->N = 0;

  A->M_padded = 0;
  A->N_padded = 0;

  A->M_block = 0;
  A->N_block = 0;
}
