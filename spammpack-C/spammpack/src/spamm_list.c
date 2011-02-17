#include "spamm_list.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/** A list. */
struct spamm_list_t
{
  /** The length of the list, i.e. the number of elements in this list. */
  unsigned int length;

  /** The elements of this list. */
  unsigned int *data;
};

/** Create a new list.
 *
 * @param length The initial size of the list.
 *
 * @return A pointer to the newly allocated and initialized list.
 */
struct spamm_list_t *
spamm_list_new (const unsigned int length)
{
  struct spamm_list_t *list;

  list = calloc(sizeof(struct spamm_list_t), 1);
  list->length = length;
  list->data = calloc(sizeof(unsigned int), length);

  return list;
}

/** Delete a list.
 *
 * @param list The list to free().
 */
void
spamm_list_delete (struct spamm_list_t **list)
{
  free((*list)->data);
  free(*list);
  *list = NULL;
}

/** Quicksort.
 */
unsigned int
spamm_list_sort_quicksort_1_partition (
    struct spamm_list_t *list,
    const unsigned int left,
    const unsigned int right,
    const unsigned int pivot,
    int (*compare) (const unsigned int, const unsigned int, void *),
    void *user_data)
{
  unsigned int i;
  unsigned int temp;
  unsigned int store_index;
  unsigned int pivot_value = list->data[pivot];

  /* Move pivot to the end. */
  temp = list->data[pivot];
  list->data[pivot] = list->data[right];
  list->data[right] = temp;

  store_index = left;

  for (i = left; i < right; i++)
  {
    if (compare(list->data[i], pivot_value, user_data) <= 0)
    {
      temp = list->data[i];
      list->data[i] = list->data[store_index];
      list->data[store_index] = temp;

      store_index++;
    }
  }

  /* Move pivot back into its final position. */
  temp = list->data[store_index];
  list->data[store_index] = list->data[right];
  list->data[right] = temp;

  return store_index;
}

/** Quicksort.
 */
void
spamm_list_sort_quicksort_1 (
    struct spamm_list_t *list,
    const unsigned int left,
    const unsigned int right,
    int (*compare) (const unsigned int, const unsigned int, void *),
    void *user_data)
{
  unsigned int pivot;
  unsigned int new_pivot;
  unsigned int length_left, length_right;

  if (right > left)
  {
    pivot = left+(right-left)/2;
    new_pivot = spamm_list_sort_quicksort_1_partition(list, left, right, pivot, compare, user_data);

    if (left+1 <= new_pivot)
    {
      length_left = new_pivot-1-left;
    }

    else { length_left = 0; }

    if (new_pivot+1 <= right)
    {
      length_right = right-new_pivot-1;
    }

    else { length_right = 0; }

    /* From http://en.wikipedia.org/wiki/Quicksort:
     *
     * Section Implementation issues:Optimizations
     *
     * To make sure at most O(log N) space is used, recurse first into the
     * smaller half of the array, and use a tail call to recurse into the
     * other.
     */
    if (length_left <= length_right)
    {
      if (length_left > 0)
      {
        spamm_list_sort_quicksort_1(list, left, new_pivot-1, compare, user_data);
      }
      if (length_right > 0)
      {
        spamm_list_sort_quicksort_1(list, new_pivot+1, right, compare, user_data);
      }
    }

    else
    {
      if (length_right > 0)
      {
        spamm_list_sort_quicksort_1(list, new_pivot+1, right, compare, user_data);
      }
      if (length_left > 0)
      {
        spamm_list_sort_quicksort_1(list, left, new_pivot-1, compare, user_data);
      }
    }
  }
}

/** Increment a somewhat signed unsigned int.
 *
 * minus_one A flag that indicates whether index is really -1.
 * index The unsigned int.
 */
void
spamm_list_sort_signed_increment (short *minus_one, unsigned int *index)
{
  if (*minus_one == 1)
  {
    *minus_one = 0;
    if (*index != 0)
    {
      printf("index should be zero\n");
      exit(1);
    }
    *index = 0;
  }

  else
  {
    (*index)++;
  }
}

/** Decrement a somewhat signed unsigned int.
 *
 * minus_one A flag that indicates whether index is really -1.
 * index The unsigned int.
 */
void
spamm_list_sort_signed_decrement (short *minus_one, unsigned int *index)
{
  if (*minus_one == 0)
  {
    if (*index == 0)
    {
      *minus_one = 1;
    }

    else
    {
      (*index)--;
    }
  }

  else
  {
    if (*index != 0)
    {
      printf("index should be zero\n");
      exit(1);
    }

    else
    {
      printf("index is already minus one, can not decrement further.\n");
      exit(1);
    }
  }
}

/** Quicksort.
 *
 * @param list The list.
 * @param left The left index of the subset to sort.
 * @param right The right index of the subset to sort.
 * @param compare The comparison function.
 * @param user_data A pointer to some data.
 */
void
spamm_list_sort_quicksort_2 (
    struct spamm_list_t *list,
    const unsigned int left,
    const unsigned int right,
    int (*compare) (const unsigned int, const unsigned int, void *),
    void *user_data)
{
  unsigned int pivot;
  unsigned int left_index = left;
  unsigned int right_index = right;
  short left_index_minus_one = 0;
  short right_index_minus_one = 0;
  unsigned int temp;

  if (right > left)
  {
    pivot = left+(right-left)/2;
    while (left_index <= pivot && pivot <= right_index)
    {
      while (compare(list->data[left_index], list->data[pivot], user_data) < 0 && left_index <= pivot)
      {
        spamm_list_sort_signed_increment(&left_index_minus_one, &left_index);
      }

      while (compare(list->data[pivot], list->data[right_index], user_data) < 0 && pivot <= right_index)
      {
        spamm_list_sort_signed_decrement(&right_index_minus_one, &right_index);
      }

      temp = list->data[left_index];
      list->data[left_index] = list->data[right_index];
      list->data[right_index] = temp;

      spamm_list_sort_signed_increment(&left_index_minus_one, &left_index);
      spamm_list_sort_signed_decrement(&right_index_minus_one, &right_index);

      if (left_index == pivot+1)
      {
        spamm_list_sort_signed_increment(&right_index_minus_one, &right_index);
        pivot = right_index;
      }

      else if ((right_index_minus_one == 1 && pivot == 0) || (right_index+1 == pivot))
      {
        spamm_list_sort_signed_decrement(&left_index_minus_one, &left_index);
        pivot = left_index;
      }
    }

    if (pivot > 0)
    {
      spamm_list_sort_quicksort_2(list, left, pivot-1, compare, user_data);
    }

    spamm_list_sort_quicksort_2(list, pivot+1, right, compare, user_data);
  }
}

/** A simple key comparison.
 *
 * @param a The first key.
 * @param b The second key.
 * @param a_norm The norm of the first key.
 * @param b_norm The norm of the second key.
 *
 * @return -1 if a < b, 0 if a == b, and 1 if a > b.
 */
int
spamm_list_compare_int (const unsigned int a, const unsigned int b, const float a_norm, const float b_norm)
{
  if (a < b)       { return -1; }
  else if (a == b) { return  0; }
  else             { return  1; }
}

/** Merge sort.
 *
 * @param list The list to sort.
 * @param scratch Scratch space, a list that has to be as long as list.
 * @param left The left most index of the sublist to sort.
 * @param right The right most index of the sublist to sort.
 * @param compare The comparison function.
 * @param user_data A pointer that will be passed to the comparison function.
 */
void
spamm_list_sort_mergesort (
    struct spamm_list_t *list,
    struct spamm_list_t *scratch,
    const unsigned int left, const unsigned int right,
    int (*compare) (const unsigned int, const unsigned int, void *),
    void *user_data)
{
  unsigned int i, i_left, i_right;
  unsigned int middle;

  if (right <= left+1) { return; }

  /* Split the list into equal halves. */
  middle = left+(right-left)/2;

  spamm_list_sort_mergesort(list, scratch, left, middle, compare, user_data);
  spamm_list_sort_mergesort(list, scratch, middle, right, compare, user_data);

  /* Merge the sorted lists. */
  for (i = left, i_left = left, i_right = middle; i < right; i++)
  {
    if (i_left < middle && i_right < right)
    {
      if (compare(list->data[i_left], list->data[i_right], user_data) <= 0)
      {
        scratch->data[i] = list->data[i_left++];
      }

      else
      {
        scratch->data[i] = list->data[i_right++];
      }
    }

    else if (i_left < middle)
    {
      scratch->data[i] = list->data[i_left++];
    }

    else
    {
      scratch->data[i] = list->data[i_right++];
    }
  }

  /* Copy the merged list back. */
  for (i = left; i < right; i++)
  {
    list->data[i] = scratch->data[i];
  }
}

/** Iterative merge sort.
 *
 * @param list The list to sort.
 * @param compare The comparison function.
 * @param user_data A pointer that will be passed to the comparison function.
 */
void
spamm_list_sort_iterative_mergesort (
    struct spamm_list_t *list,
    float *node_norm,
    spamm_list_compare_function compare)
{
  unsigned int i, j, j_next, i_left, i_right;
  unsigned int sub_current, sub_next;
  unsigned int sub_length;
  unsigned int *sublist;
  unsigned int *scratch;
  float *scratch_norm;

  /* The list is trivially already sorted. */
  if (list->length <= 1) { return; }

  /* Create index array for sublists. This array is length+1 since we
   * terminate the array by a value of length. */
  sublist = (unsigned int*) malloc(sizeof(unsigned int)*2*(list->length+1));

  /* Allocate scratch space. */
  scratch = (unsigned int*) calloc(sizeof(unsigned int), list->length);
  scratch_norm = (float*) calloc(sizeof(float), list->length);

  /* Break the original list into at most N pieces, i.e. single element
   * sublists. If adajacent list elements are already in the right order, we
   * put them into the same sublist. */
  sublist[0] = 0;
  for (i = 1, j = 1; i < list->length; i++)
  {
    if (compare(list->data[i-1], list->data[i], node_norm[i-1], node_norm[i]) > 0)
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
  sub_next = list->length+1;
  while (sublist[sub_current+1] < list->length)
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
          if (compare(list->data[i_left], list->data[i_right], node_norm[i_left], node_norm[i_right]) <= 0)
          {
            scratch[i] = list->data[i_left];
            scratch_norm[i] = node_norm[i_left];
            i_left++;
          }

          else
          {
            scratch[i] = list->data[i_right];
            scratch_norm[i] = node_norm[i_right];
            i_right++;
          }
        }

        else if (i_left < sublist[sub_current+j+1])
        {
          scratch[i] = list->data[i_left];
          scratch_norm[i] = node_norm[i_left];
          i_left++;
        }

        else
        {
          scratch[i] = list->data[i_right];
          scratch_norm[i] = node_norm[i_right];
          i_right++;
        }
      }

      /* Copy the merged list back. */
      for (i = sublist[sub_current+j]; i < sublist[sub_current+j+2]; i++)
      {
        list->data[i] = scratch[i];
        node_norm[i] = scratch_norm[i];
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
    if (sub_current == 0) { sub_current = list->length+1; }
    else                  { sub_current = 0; }
    if (sub_next == 0) { sub_next = list->length+1; }
    else               { sub_next = 0; }
  }

  /* Free memory. */
  free(sublist);
  free(scratch);
  free(scratch_norm);
}

/** Sort a list.
 *
 * @param list The list to sort.
 * @param node_norm The node norms.
 * @param compare A function that compares 2 elements a and b of the list, and
 * returns -1 in case a < b, 0 if a == b, and 1 if a > b.
 */
void
spamm_list_sort (
    struct spamm_list_t *list,
    float *node_norm,
    spamm_list_compare_function compare)
{
  struct spamm_list_t *scratch;

#ifdef SPAMM_SORT_QUICKSORT
  spamm_list_sort_quicksort_1(list, 0, list->length-1, compare, user_data);
#endif

#ifdef SPAMM_SORT_MERGESORT
  scratch = spamm_list_new(list->length);
  spamm_list_sort_mergesort(list, scratch, 0, list->length, compare, user_data);
  spamm_list_delete(&scratch);
#endif

  spamm_list_sort_iterative_mergesort(list, node_norm, compare);
}

/** Return the length of a list.
 *
 * @param list The list.
 *
 * @return The number of elements stored in the list.
 */
unsigned int
spamm_list_length (struct spamm_list_t *list)
{
  return list->length;
}

/** Get an element from a list.
 *
 * @param list The list.
 * @param i The index of the entry to get.
 *
 * @return The entry.
 */
unsigned int
spamm_list_get (struct spamm_list_t *list, const unsigned int i)
{
  assert(list != NULL);

  if (i >= list->length)
  {
    printf("index out of bounds\n");
    exit(1);
  }

  return list->data[i];
}

/** Set an element in the list.
 *
 * @param list The list.
 * @param i The index of the entry to set.
 * @param Ai The value of the entry.
 */
void
spamm_list_set (struct spamm_list_t *list, const unsigned int i, const unsigned int Ai)
{
  assert(list != NULL);

  if (i >= list->length)
  {
    printf("index out of bounds\n");
    exit(1);
  }

  list->data[i] = Ai;
}

/** Get the pointer to the data.
 *
 * @param list The list.
 *
 * @return A pointer unsigned int* to the elements in this list.
 */
unsigned int*
spamm_list_get_data (struct spamm_list_t *list)
{
  return list->data;
}
