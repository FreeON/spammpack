#include "spamm.h"

/** Compare two spamm_multiply_stream_element_t.
 *
 * @param element1 The first stream element.
 * @param element2 The second stream element.
 *
 * @return element1->C_index < element2->C_index: -1, element1->C_index ==
 *         element2->C_index: 0, element1->C_index > element2->C_index: +1.
 */
int
spamm_compare_multiply_stream_element (const void *element1, const void *element2)
{
  const struct spamm_multiply_stream_element_t *ele1 = element1;
  const struct spamm_multiply_stream_element_t *ele2 = element2;

  if (ele1->C_index < ele2->C_index) { return -1; }
  else if (ele1->C_index == ele2->C_index) { return 0; }
  else { return 1; }
}
