/** @file */

#include "spamm_sort.h"

#include <stdlib.h>

/** Iterative merge sort.
 *
 * @param length The number of elements in the list.
 * @param list The list to sort.
 * @param mask The mask to apply to the values in the list before sorting.
 */
void
#ifdef SPAMM_SORT_MASKED
SPAMM_FUNCTION(spamm_sort_masked, SPAMM_FUNC_TYPE) (const unsigned int length,
    SPAMM_SORT_TYPE *list,
    const unsigned int mask)
#else
SPAMM_FUNCTION(spamm_sort, SPAMM_FUNC_TYPE) (const unsigned int length,
    SPAMM_SORT_TYPE *list)
#endif
{
  unsigned int i, j, j_next, i_left, i_right;
  unsigned int sub_current, sub_next;
  unsigned int sub_length;
  unsigned int *sublist;
  SPAMM_SORT_TYPE *scratch;

  /* The list is trivially already sorted. */
  if (length <= 1) { return; }

  /* Create index array for sublists. This array is length+1 since we
   * terminate the array by a value of length. */
  sublist = (unsigned int*) malloc(sizeof(unsigned int)*2*(length+1));

  /* Allocate scratch space. */
  scratch = (SPAMM_SORT_TYPE*) calloc(sizeof(SPAMM_SORT_TYPE), length);

  /* Break the original list into at most N pieces, i.e. single element
   * sublists. If adajacent list elements are already in the right order, we
   * put them into the same sublist. */
  sublist[0] = 0;
  for (i = 1, j = 1; i < length; i++)
  {
#ifdef SPAMM_SORT_MASKED
    if (SPAMM_SORT_COMPARE(list[i-1], list[i], mask) > 0)
#else
    if (list[i-1] > list[i])
#endif
    {
      /* The 2 elements are in incorrect order. Start a new sublist. */
      sublist[j++] = i;
    }
  }
  sublist[j++] = i;
  sub_length = j;

  /* Loop over the list, merging neighboring sublists until everthying is
   * sorted. */
  sub_current = 0;
  sub_next = length+1;
  while (sublist[sub_current+1] < length)
  {
    for (j = 0, j_next = 0; j < sub_length-2; j += 2)
    {
      /* Merge 2 adjacent sublists. */
      for (i = sublist[sub_current+j], i_left = sublist[sub_current+j], i_right = sublist[sub_current+j+1];
          i < sublist[sub_current+j+2];
          i++)
      {
        if (i_left < sublist[sub_current+j+1] && i_right < sublist[sub_current+j+2])
        {
#ifdef SPAMM_SORT_MASKED
          if (SPAMM_SORT_COMPARE(list[i_left], list[i_right], mask) <= 0)
#else
          if (list[i_left] <= list[i_right])
#endif
          {
            scratch[i] = list[i_left];
            i_left++;
          }

          else
          {
            scratch[i] = list[i_right];
            i_right++;
          }
        }

        else if (i_left < sublist[sub_current+j+1])
        {
          scratch[i] = list[i_left];
          i_left++;
        }

        else
        {
          scratch[i] = list[i_right];
          i_right++;
        }
      }

      /* Copy the merged list back. */
      for (i = sublist[sub_current+j]; i < sublist[sub_current+j+2]; i++)
      {
        list[i] = scratch[i];
      }

      /* Remove division between the sublists just merged. */
      sublist[sub_next+j_next] = sublist[sub_current+j];
      sublist[sub_next+j_next+1] = sublist[sub_current+j+2];
      j_next++;
    }

    /* Add remaining sublist divisions. */
    while (j < sub_length)
    {
      sublist[sub_next+j_next] = sublist[sub_current+j];
      j++;
      j_next++;
    }
    sub_length = j_next;

    /* Switch sublists. */
    if (sub_current == 0) { sub_current = length+1; }
    else                  { sub_current = 0; }
    if (sub_next == 0) { sub_next = length+1; }
    else               { sub_next = 0; }
  }

  /* Free memory. */
  free(sublist);
  free(scratch);
}
