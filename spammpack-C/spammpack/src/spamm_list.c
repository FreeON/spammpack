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
spamm_list_sort_quicksort_1_partition (struct spamm_list_t *list,
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
spamm_list_sort_quicksort_1 (struct spamm_list_t *list,
    const unsigned int left, const unsigned int right,
    int (*compare) (const unsigned int, const unsigned int, void *),
    void *user_data)
{
  unsigned int pivot;
  unsigned int new_pivot;

  if (right > left)
  {
    pivot = left+(right-left)/2;
    new_pivot = spamm_list_sort_quicksort_1_partition(list, left, right, pivot, compare, user_data);

    if (left+1 <= new_pivot)
    {
      spamm_list_sort_quicksort_1(list, left, new_pivot-1, compare, user_data);
    }

    if (new_pivot+1 <= right)
    {
      spamm_list_sort_quicksort_1(list, new_pivot+1, right, compare, user_data);
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
spamm_list_sort_quicksort_2 (struct spamm_list_t *list,
    const unsigned int left, const unsigned int right,
    int (*compare) (const unsigned int, const unsigned int, void *),
    void *user_data)
{
  unsigned int pivot;
  unsigned int left_index = left;
  unsigned int right_index = right;

  unsigned int temp;

  if (right > left)
  {
    pivot = (left+right)/2;
    while (left_index <= pivot && pivot <= right_index)
    {
      while (compare(list->data[left_index], list->data[pivot], user_data) < 0 && left_index <= pivot)
      {
        left_index++;
      }

      while (compare(list->data[pivot], list->data[right_index], user_data) < 0 && pivot <= right_index)
      {
        right_index--;
      }

      temp = list->data[left_index];
      list->data[left_index] = list->data[right_index];
      list->data[right_index] = temp;

      if (left_index < right)
      {
        left_index++;
      }

      if (right_index > left)
      {
        right_index--;
      }

      if (left_index == pivot+1)
      {
        pivot = right_index;
        right_index++;
      }

      else if (right_index+1 == pivot)
      {
        pivot = left_index;
        left_index--;
      }
    }

    if (left+1 <= pivot)
    {
      spamm_list_sort_quicksort_2(list, left, pivot-1, compare, user_data);
    }

    if (pivot+1 <= right)
    {
      spamm_list_sort_quicksort_2(list, pivot+1, right, compare, user_data);
    }
  }
}

/** A simple key comparison.
 *
 * @param a The first key.
 * @param b The second key.
 * @param user_data A pointer to some data to pass to the comparison function.
 *
 * @return -1 if a < b, 0 if a == b, and 1 if a > b.
 */
int
spamm_list_compare_int (const unsigned int a, const unsigned int b, void *user_data)
{
  if (a < b)       { return -1; }
  else if (a == b) { return  0; }
  else             { return  1; }
}

/** Sort a list.
 *
 * @param list The list to sort.
 * @param compare A function that compares 2 elements a and b of the list, and
 * returns -1 in case a < b, 0 if a == b, and 1 if a > b.
 * @param user_data A pointer that is given to the compare function and can be
 * used to pass any information to the compare function.
 */
void
spamm_list_sort (struct spamm_list_t *list,
    int (*compare) (const unsigned int, const unsigned int, void *),
    void *user_data)
{
  spamm_list_sort_quicksort_1(list, 0, list->length-1, compare, user_data);
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
