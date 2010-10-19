#include "spamm.h"
#include <stdlib.h>

/** Delete a struct spamm_multiply_stream_element_t.
 *
 * @param data The pointer to the spamm_multiply_stream_element_t to delete.
 */
void
spamm_delete_multiply_stream_element (void *data)
{
  struct spamm_multiply_stream_element_t *element = data;

  if (element->C_node->block_dense != element->C_block_dense)
  {
    /* Assume that this block was allocated separately as redundant C sub
     * matrix block.
     */
    free(element->C_block_dense);
  }
  free(element);
}
