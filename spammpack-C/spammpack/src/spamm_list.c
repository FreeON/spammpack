#include "spamm_list.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/** 01010101010101010101010101010101 = 0x55555555 */
#define MASK_2D_J  0x55555555

/** 10101010101010101010101010101010 = 0xaaaaaaaa */
#define MASK_2D_I  0xaaaaaaaa

/** A list of matrix indices and their norms. */
struct spamm_list_t
{
  /** The length of the list, i.e. the number of elements in this list. */
  unsigned int length;

  /** The matrix indices. */
  unsigned int *index;

  /** The norms. */
  float *norm;
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
  list->index = (unsigned int*) calloc(sizeof(unsigned int), length);
  list->norm = (float*) calloc(sizeof(float), length);

  return list;
}

/** Delete a list.
 *
 * @param list The list to free().
 */
void
spamm_list_delete (struct spamm_list_t **list)
{
  free((*list)->index);
  free((*list)->norm);
  free(*list);
  *list = NULL;
}

/** A simple key comparison.
 *
 * @param a The first key.
 * @param b The second key.
 *
 * @return -1 if a < b, 0 if a == b, and 1 if a > b.
 */
int
spamm_list_compare_int (const unsigned int a, const unsigned int b)
{
  if (a < b)       { return -1; }
  else if (a == b) { return  0; }
  else             { return  1; }
}

/** Compare 2 2D indices by their row index.
 *
 * @param a The first index.
 * @param b The second index.
 *
 * @return if a is before b, return -1, if a is after b, return +1, and if a
 * and b are equivalent, return 0.
 */
int
spamm_list_compare_index_row (const unsigned int a, const unsigned int b)
{
  unsigned int a_masked = a & MASK_2D_I;
  unsigned int b_masked = b & MASK_2D_I;

  if (a_masked < b_masked)       { return -1; }
  else if (a_masked > b_masked)  { return  1; }
  else                           { return 0; }
}

/** Compare 2 2D indices by their column index.
 *
 * @param a The first index.
 * @param b The second index.
 *
 * @return if a is before b, return -1, if a is after b, return +1, and if a
 * and b are equivalent, return 0.
 */
int
spamm_list_compare_index_column (const unsigned int a, const unsigned int b)
{
  unsigned int a_masked = a & MASK_2D_J;
  unsigned int b_masked = b & MASK_2D_J;

  if (a_masked < b_masked)       { return -1; }
  else if (a_masked > b_masked)  { return  1; }
  else                           { return 0; }
}

/** Iterative merge sort.
 *
 * @param list The list to sort.
 * @param compare The comparison function.
 * @param user_data A pointer that will be passed to the comparison function.
 */
void
spamm_list_sort_index_iterative_mergesort (
    struct spamm_list_t *list,
    spamm_list_compare_index_function compare)
{
  unsigned int i, j, j_next, i_left, i_right;
  unsigned int sub_current, sub_next;
  unsigned int sub_length;
  unsigned int *sublist;
  unsigned int *scratch_index;
  float *scratch_norm;

  /* The list is trivially already sorted. */
  if (list->length <= 1) { return; }

  /* Create index array for sublists. This array is length+1 since we
   * terminate the array by a value of length. */
  sublist = (unsigned int*) malloc(sizeof(unsigned int)*2*(list->length+1));

  /* Allocate scratch space. */
  scratch_index = (unsigned int*) calloc(sizeof(unsigned int), list->length);
  scratch_norm = (float*) calloc(sizeof(float), list->length);

  /* Break the original list into at most N pieces, i.e. single element
   * sublists. If adajacent list elements are already in the right order, we
   * put them into the same sublist. */
  sublist[0] = 0;
  for (i = 1, j = 1; i < list->length; i++)
  {
    if (compare(list->index[i-1], list->index[i]) > 0)
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
          if (compare(list->index[i_left], list->index[i_right]) <= 0)
          {
            scratch_index[i] = list->index[i_left];
            scratch_norm[i] = list->norm[i_left];
            i_left++;
          }

          else
          {
            scratch_index[i] = list->index[i_right];
            scratch_norm[i] = list->norm[i_right];
            i_right++;
          }
        }

        else if (i_left < sublist[sub_current+j+1])
        {
          scratch_index[i] = list->index[i_left];
          scratch_norm[i] = list->norm[i_left];
          i_left++;
        }

        else
        {
          scratch_index[i] = list->index[i_right];
          scratch_norm[i] = list->norm[i_right];
          i_right++;
        }
      }

      /* Copy the merged list back. */
      for (i = sublist[sub_current+j]; i < sublist[sub_current+j+2]; i++)
      {
        list->index[i] = scratch_index[i];
        list->norm[i] = scratch_norm[i];
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
  free(scratch_index);
  free(scratch_norm);
}

/** Iterative merge sort.
 *
 * @param list The list to sort.
 * @param compare The comparison function.
 * @param user_data A pointer that will be passed to the comparison function.
 */
void
spamm_list_sort_norm_iterative_mergesort (
    struct spamm_list_t *list,
    const unsigned int left,
    const unsigned int right)
{
  unsigned int i, j, j_next, i_left, i_right;
  unsigned int length;
  unsigned int sub_current, sub_next;
  unsigned int sub_length;
  unsigned int *sublist;
  unsigned int *scratch_index;
  float *scratch_norm;

  /* What is the length of this sublist? */
  if (right >= left) { length = right-left; }
  else { return; }

  /* The list is trivially already sorted. */
  if (length <= 1) { return; }

  /* Create index array for sublists. This array is length+1 since we
   * terminate the array by a value of length. */
  sublist = (unsigned int*) malloc(sizeof(unsigned int)*2*(length+1));

  /* Allocate scratch space. */
  scratch_index = (unsigned int*) malloc(sizeof(unsigned int)*length);
  scratch_norm = (float*) malloc(sizeof(float)*length);

  /* Break the original list into at most N pieces, i.e. single element
   * sublists. If adjacent list elements are already in the right order, we
   * put them into the same sublist. */
  sublist[0] = left;
  for (i = left+1, j = 1; i < right; i++)
  {
    if (list->norm[i-1] < list->norm[i])
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
  while (sublist[sub_current+1] < right)
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
          if (list->norm[i_left] >= list->norm[i_right])
          {
            /* The left element is in the right place. */
            scratch_index[i-left] = list->index[i_left];
            scratch_norm[i-left] = list->norm[i_left];
            i_left++;
          }

          else
          {
            /* The right element is in the right place. */
            scratch_index[i-left] = list->index[i_right];
            scratch_norm[i-left] = list->norm[i_right];
            i_right++;
          }
        }

        /* Merge the left over elements. */
        else if (i_left < sublist[sub_current+j+1])
        {
          scratch_index[i-left] = list->index[i_left];
          scratch_norm[i-left] = list->norm[i_left];
          i_left++;
        }

        else
        {
          scratch_index[i-left] = list->index[i_right];
          scratch_norm[i-left] = list->norm[i_right];
          i_right++;
        }
      }

      /* Copy the merged list back. */
      for (i = sublist[sub_current+j]; i < sublist[sub_current+j+2]; i++)
      {
        list->index[i] = scratch_index[i-left];
        list->norm[i] = scratch_norm[i-left];
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
  free(scratch_index);
  free(scratch_norm);
}

/** Sort a list by its matrix index.
 *
 * @param list The list to sort.
 * @param compare A function that compares 2 elements a and b of the list, and
 * returns -1 in case a < b, 0 if a == b, and 1 if a > b.
 */
void
spamm_list_sort_index (struct spamm_list_t *list, spamm_list_compare_index_function compare)
{
  spamm_list_sort_index_iterative_mergesort(list, compare);
}

/** Sort a list by its matrix norm.
 *
 * @param list The list to sort.
 * @param left The left most index of the sub-list.
 * @param right The right most index of the sub-list.
 */
void
spamm_list_sort_norm (struct spamm_list_t *list, const unsigned int left, const unsigned int right)
{
  spamm_list_sort_norm_iterative_mergesort(list, left, right);
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

/** Get an index from a list.
 *
 * @param list The list.
 * @param i The index of the entry to get.
 *
 * @return The entry (the matrix index).
 */
unsigned int
spamm_list_get_index (struct spamm_list_t *list, const unsigned int i)
{
  assert(list != NULL);

  if (i >= list->length)
  {
    printf("[list get index] index out of bounds: (i = %u) >= (length = %u)\n", i, list->length);
    exit(1);
  }

  return list->index[i];
}

/** Get a norm from a list.
 *
 * @param list The list.
 * @param i The index of the entry to get.
 *
 * @return The entry (the matrix norm).
 */
float
spamm_list_get_norm (struct spamm_list_t *list, const unsigned int i)
{
  assert(list != NULL);

  if (i >= list->length)
  {
    printf("[list get norm] index out of bounds: (i = %u) >= (length = %u)\n", i, list->length);
    exit(1);
  }

  return list->norm[i];
}

/** Set an element in the list.
 *
 * @param list The list.
 * @param i The index of the entry to set.
 * @param index The index of the entry.
 * @param norm The norm of the entry.
 */
void
spamm_list_set (struct spamm_list_t *list, const unsigned int i, const unsigned int index, const float norm)
{
  assert(list != NULL);

  if (i >= list->length)
  {
    printf("[list set] index out of bounds: (i = %u) >= (length = %u)\n", i, list->length);
    exit(1);
  }

  list->index[i] = index;
  list->norm[i] = norm;
}

/** Print a list.
 *
 * @param list The list.
 */
void
spamm_list_print (const struct spamm_list_t *list)
{
  unsigned int i;

  printf("# list with %u entries\n", list->length);
  printf("# i index column row norm\n");
  for (i = 0; i < list->length; i++)
  {
    printf("%i %i %i %i %e\n", i, list->index[i], list->index[i] & MASK_2D_J, list->index[i] & MASK_2D_I, list->norm[i]);
  }
}
