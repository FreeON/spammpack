#include "spamm.h"
#include <stdio.h>
#include <stdlib.h>

/** Print a multiply stream.
 *
 * @param stream The multiply stream.
 */
void
spamm_print_multiply_stream (const struct spamm_ll_t *stream)
{
  unsigned int i;
  struct spamm_ll_iterator_t *iterator;
  struct spamm_ll_node_t *node;
  struct spamm_multiply_stream_element_t *element;

  iterator = spamm_ll_iterator_new(stream);
  printf("multiply stream (%u elements):", stream->number_elements);
  for (i = 0, node = spamm_ll_iterator_first(iterator); node != NULL; node = spamm_ll_iterator_next(iterator))
  {
    element = node->data;
    printf(" (%u) index = %u:%u:%u, block = %f", i++, element->A_index, element->B_index, element->C_index, element->C_block_dense[0]);
  }
  printf("\n");

  spamm_ll_iterator_delete(&iterator);
}
