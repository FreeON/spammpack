/** @file */

#include "spamm.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/** Initialize a new linked list.
 *
 * @param list The linked list to initialize.
 */
void
spamm_ll_new (struct spamm_multiply_stream_t *list)
{
  assert(list != NULL);

  list->number_elements = 0;
  list->first = NULL;
  list->last = NULL;
}

/** Delete a linked list.
 *
 * @param list The linked list to delete.
 */
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

/** Initialize a new node of a linked list.
 *
 * @param node The linked list node to initialize.
 */
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

  (*node)->alpha = 0;
  (*node)->beta = 0;

  (*node)->A_node = NULL;
  (*node)->B_node = NULL;
  (*node)->C_node = NULL;
}

/** Append a node to a linked list.
 *
 * @param alpha The value of alpha of spamm_node_t.
 * @param beta The value of beta of spamm_node_t.
 * @param A_index The linear index of A_node in spamm_node_t.
 * @param A_node The spamm_node_t A_node.
 * @param B_index The linear index of B_node in spamm_node_t.
 * @param B_node The spamm_node_t B_node.
 * @param C_index The linear index of C_node in spamm_node_t.
 * @param C_node The spamm_node_t C_node.
 * @param list The linked list.
 */
void
spamm_ll_append (const float_t alpha, const float_t beta,
    const unsigned int A_index, struct spamm_node_t *A_node,
    const unsigned int B_index, struct spamm_node_t *B_node,
    const unsigned int C_index, struct spamm_node_t *C_node,
    struct spamm_multiply_stream_t *list)
{
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

  new_node->alpha = alpha;
  new_node->beta = beta;

  /* Append new node to list. */
  new_node->previous = list->last;
  if (list->first == NULL) { list->first = new_node; }
  if (list->last != NULL) { list->last->next = new_node; }
  list->last = new_node;

  /* Increment element counter. */
  list->number_elements++;
}

/** Get an element from a linked list.
 *
 * @param i Index of element to get. Counting starts with 0.
 * @param list The linked list.
 *
 * @return The element of the linked list. A return value of NULL means that
 *         the element was not found in the list.
 */
struct spamm_multiply_stream_node_t *
spamm_ll_get (const unsigned int i, const struct spamm_multiply_stream_t *list)
{
  int j;
  struct spamm_multiply_stream_node_t *node;

  assert(list != NULL);

  if (i >= list->number_elements)
  {
    return NULL;
  }

  j = 0;
  for (node = list->first; node != NULL; node = node->next)
  {
    if (i == j) { break; }
    j++;
  }

  if (j < list->number_elements && node == NULL)
  {
    spamm_log("bug? i = %i, j = %i\n", __FILE__, __LINE__, i, j);
    spamm_ll_print(list);
    exit(1);
  }

  return node;
}

/** Swap to elements in a linked list.
 *
 * @param node1 The first node of the list.
 * @param node2 The second node of the list.
 * @param list The linked list.
 */
void
spamm_ll_swap (struct spamm_multiply_stream_node_t **node1,
    struct spamm_multiply_stream_node_t **node2,
    struct spamm_multiply_stream_t *list)
{
  struct spamm_multiply_stream_node_t *node1_previous = (*node1)->previous;
  struct spamm_multiply_stream_node_t *node2_previous = (*node2)->previous;
  struct spamm_multiply_stream_node_t *node1_next = (*node1)->next;
  struct spamm_multiply_stream_node_t *node2_next = (*node2)->next;

  struct spamm_multiply_stream_node_t *temp;

  assert(list != NULL);
  assert(*node1 != NULL);
  assert(*node2 != NULL);

  /* Swap out first and last link. */
  if      (list->first == *node1) { list->first = *node2; }
  else if (list->first == *node2) { list->first = *node1; }
  if      (list->last == *node1)  { list->last = *node2; }
  else if (list->last == *node2)  { list->last = *node1; }

  /* Connect neighbors of node1 and node2. */
  if (node1_previous != NULL) { node1_previous->next = *node2; }
  if (node2_previous != NULL) { node2_previous->next = *node1; }
  if (node1_next != NULL) { node1_next->previous = *node2; }
  if (node2_next != NULL) { node2_next->previous = *node1; }

  /* Connect to neighbors. */
  temp = (*node1)->previous; (*node1)->previous = (*node2)->previous; (*node2)->previous = temp;
  temp = (*node1)->next; (*node1)->next = (*node2)->next; (*node2)->next = temp;

  /* Swap nodes. */
  temp = *node1; *node1 = *node2; *node2 = temp;
}

/** Sort a linked list.
 *
 * @param list The linked list.
 */
void
spamm_ll_sort (struct spamm_multiply_stream_t *list)
{
  /* We optimize locality by sorting the stream so as to minimize the change
   * of any of the indices. */

  struct spamm_multiply_stream_node_t *node1, *node2;

  assert(list != NULL);

  for (node1 = list->first; node1 != list->last && node1 != NULL; node1 = node1->next) {
    for (node2 = node1->next; node2 != NULL; node2 = node2->next)
    {
      if ((node1->A_index+node1->B_index+node1->C_index) > (node2->A_index+node2->B_index+node2->C_index))
      {
        spamm_ll_swap(&node1, &node2, list);
      }
    }
  }
}

/** Print debugging information of a linked list node.
 *
 * @param name A string that gets printed as the name of the node.
 * @param node The linked list node to print.
 */
void
spamm_ll_print_node_debug (const char *name, const struct spamm_multiply_stream_node_t *node)
{
  assert(node != NULL);

  printf("%s (%p): previous = %p, next = %p, (%u,%u,%u), alpha = %f, beta = %f\n",
      name, node, node->previous, node->next, node->A_index, node->B_index, node->C_index,
      node->alpha, node->beta);
}

/** Print information of a linked list node.
 *
 * @param node The linked list node.
 */
void
spamm_ll_print_node (const struct spamm_multiply_stream_node_t *node)
{
  assert(node != NULL);

  printf("(%u,%u,%u)\n", node->A_index, node->B_index, node->C_index);
}

/** Print a linked list.
 *
 * @param list The linked list to print.
 */
void
spamm_ll_print (const struct spamm_multiply_stream_t *list)
{
  unsigned int i;

  assert(list != NULL);

  struct spamm_multiply_stream_node_t *node;

  printf("multiply stream: [ ");
  i = 0;
  for (node = list->first; node != NULL; node = node->next)
  {
    printf("%u:(%u,%u,%u) ", ++i, node->A_index, node->B_index, node->C_index);
  }
  printf("]\n");
}

/** Print a linked list in matlab format.
 *
 * @param list The linked list.
 */
void
spamm_ll_print_matlab (const struct spamm_multiply_stream_t *list)
{
  assert(list != NULL);

  struct spamm_multiply_stream_node_t *node;

  printf("MStream = [\n");
  for (node = list->first; node != NULL; node = node->next)
  {
    printf("[%u,%u,%u]\n", node->A_index, node->B_index, node->C_index);
  }
  printf("];\n");
}
