#include "spamm.h"

#include <stdlib.h>

/** Delete a block of type spamm_data_t.
 *
 * @param data The block of data to free.
 */
void
spamm_delete_block (struct spamm_data_t **data)
{
  free(*data);
  *data = NULL;
}