/** @file */

#include <stdlib.h>

/** Compare 2 2D indices with a mask.
 *
 * @param a The first index.
 * @param b The second index.
 * @param mask The mask.
 *
 * @return if a is before b, return -1, if a is after b, return +1, and if a
 * and b are equivalent, return 0.
 */
int
spamm_compare_masked (const unsigned int a,
    const unsigned int b,
    const unsigned int mask)
{
  unsigned int a_masked = a & mask;
  unsigned int b_masked = b & mask;

  if (a_masked < b_masked)      { return -1; }
  else if (a_masked > b_masked) { return  1; }
  else                          { return  0; }
}

/** Iterative merge sort.
 *
 * @param length The number of elements in the list.
 * @param list The list to sort.
 * @param mask The mask to apply to the values in the list before sorting.
 */
void
SPAMM_FUNCTION(spamm_sort_masked, SPAMM_FUNC_TYPE) (const unsigned int length,
    SPAMM_SORT_TYPE *list,
    const unsigned int mask)
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
  scratch = (unsigned int*) calloc(sizeof(unsigned int), length);

  /* Break the original list into at most N pieces, i.e. single element
   * sublists. If adajacent list elements are already in the right order, we
   * put them into the same sublist. */
  sublist[0] = 0;
  for (i = 1, j = 1; i < length; i++)
  {
    if (spamm_compare_masked(list[i-1], list[i], mask) > 0)
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
          if (spamm_compare_masked(list[i_left], list[i_right], mask) <= 0)
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
