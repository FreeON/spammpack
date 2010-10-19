#include "spamm.h"
#include <stdlib.h>

/** Mask index i or j out of linear index.
 *
 * The parameter index is interpreted as a bit field representing the
 * hierarchical tree structure. The masking operation removes either the i or
 * the j component and returns the linear index with the remaining component.
 * As an example consider the linear index:
 *
 * index = 0100110100
 *
 * An operation spamm_mask(index, mask_i) removes the j-component. The result
 * is:
 *
 * \code
 * spamm_mask(index = 0100110100, mask_i) = 0000100000
 * \endcode
 *
 * Applying the other mask yields:
 *
 * \code
 * spamm_mask(index = 0100110100, mask_j) = 0100010100
 * \endcode
 *
 * and the original index is recovered by OR'ing the two masked indices:
 *
 * \code
 *   0000100000
 * | 0100010100
 * ------------
 *   0100110100
 * \endcode
 *
 * @param index The linear index.
 * @param width The number of bits to consider in the masking operation.
 * @param mask The mask.
 *
 * @return The masked linear index.
 */
unsigned int
spamm_mask (const unsigned int index, const unsigned int width,
    const enum spamm_linear_mask_t mask)
{
  int i;
  unsigned int bitmask = 0;
  unsigned int result = 0;

  switch (mask)
  {
    case i_mask:
      bitmask = 2;
      break;

    case j_mask:
      bitmask = 1;
      break;

    default:
      LOG2_FATAL("unknown mask\n");
      exit(1);
      break;
  }

  for (i = 0; i < width; i += 2)
  {
    if ((index & bitmask) != 0)
    {
      result |= bitmask;
    }

    bitmask <<= 2;
  }

  return result;
}
