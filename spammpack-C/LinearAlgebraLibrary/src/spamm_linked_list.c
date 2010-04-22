#include "spamm.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

void
spamm_ll_new (struct spamm_multiply_stream_t *list)
{
  assert(list != NULL);

  list->first = NULL;
  list->last = NULL;
}

void
spamm_ll_delete (struct spamm_multiply_stream_t *list)
{
  struct spamm_multiply_stream_node_t *node1, *node2;

  assert(list != NULL);

  for (node1 = list->first; node1 != NULL; )
  {
    node2 = node1;
    node1 = node1->next;

    free(node2);
  }
  spamm_ll_new(list);
}

void
spamm_ll_new_node (struct spamm_multiply_stream_node_t **node)
{
  assert(node != NULL);

  *node = (struct spamm_multiply_stream_node_t*) malloc(sizeof(struct spamm_multiply_stream_node_t));

  (*node)->previous = NULL;
  (*node)->next = NULL;

  (*node)->A_index = 0;
  (*node)->B_index = 0;
  (*node)->C_index = 0;

  (*node)->A_node = NULL;
  (*node)->B_node = NULL;
  (*node)->C_node = NULL;
}

void
spamm_ll_append (const unsigned int A_index, const struct spamm_node_t *A_node,
    const unsigned int B_index, const struct spamm_node_t *B_node,
    const unsigned int C_index, const struct spamm_node_t *C_node,
    struct spamm_multiply_stream_t *list)
{
  struct spamm_multiply_stream_node_t *node;
  struct spamm_multiply_stream_node_t *new_node;

  assert(list != NULL);
  assert(A_node != NULL);
  assert(B_node != NULL);
  assert(C_node != NULL);

  spamm_ll_new_node(&new_node);

  new_node->A_index = A_index;
  new_node->B_index = B_index;
  new_node->C_index = C_index;
  new_node->A_node  = A_node;
  new_node->B_node  = B_node;
  new_node->C_node  = C_node;

  /* Append new node to list. */
  new_node->previous = list->last;
  if (list->first == NULL) { list->first = new_node; }
  if (list->last != NULL) { list->last->next = new_node; }
  list->last = new_node;
}

void
spamm_ll_sort (struct spamm_multiply_stream_t *list)
{
}

void
spamm_ll_print (struct spamm_multiply_stream_t *list)
{
  assert(list != NULL);

  struct spamm_multiply_stream_node_t *node;

  printf("list: [ ");
  for (node = list->first; node != NULL; node = node->next)
  {
    printf("(%u,%u,%u) ", node->A_index, node->B_index, node->C_index);
  }
  printf("]\n");
}
