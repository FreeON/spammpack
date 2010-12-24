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
    int (*compare) (const void *, const void *, void *),
    void *user_data)
{
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
