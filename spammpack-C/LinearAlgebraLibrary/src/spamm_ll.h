/** @file */

#if ! defined(__SPAMM_LL_H)

/** Define in case spamm_ll.h has been included. */
#define __SPAMM_LL_H 1

/** \page page_ll The SpAMM linked list.
 *
 * Some functions for dealing with linked lists. Details can be found in
 * spamm_ll.h.
 */

/** The linked list.
 */
struct spamm_ll_t
{
  /** Number of elements. */
  unsigned int number_elements;

  /** Links to the first node in list. */
  struct spamm_ll_node_t *first;

  /** Links to the last node in list. */
  struct spamm_ll_node_t *last;
};

/** A node in a linked list.
 */
struct spamm_ll_node_t
{
  /** Link to the previous element. */
  struct spamm_ll_node_t *previous;

  /** Link to the next element. */
  struct spamm_ll_node_t *next;

  /** Pointer to node data. */
  void *data;
};

/** An iterator object. */
struct spamm_ll_iterator_t
{
  /** The list this iterator is applied to. */
  struct spamm_ll_t *list;

  /** The last node returned. */
  struct spamm_ll_node_t *last_node;
};

void
spamm_ll_append (void *data, struct spamm_ll_t *list);

void
spamm_ll_delete (struct spamm_ll_t *list);

void *
spamm_ll_get (const unsigned int i, const struct spamm_ll_t *list);

void
spamm_ll_iterator_delete (struct spamm_ll_iterator_t **iterator);

struct spamm_ll_node_t *
spamm_ll_iterator_first (struct spamm_ll_iterator_t *iterator);

struct spamm_ll_iterator_t *
spamm_ll_iterator_new (struct spamm_ll_t *list);

struct spamm_ll_node_t *
spamm_ll_iterator_next (struct spamm_ll_iterator_t *iterator);

struct spamm_ll_t *
spamm_ll_new ();

struct spamm_ll_node_t *
spamm_ll_new_node ();

void
spamm_ll_print (char *(*data_to_string) (const void *data), const struct spamm_ll_t *list);

void
spamm_ll_print_node (char *(*data_to_string) (const void *data), const struct spamm_ll_node_t *node);

void
spamm_ll_sort (int (*compare) (const void *data1, const void *data2), struct spamm_ll_t *list);

void
spamm_ll_swap (struct spamm_ll_node_t **node1, struct spamm_ll_node_t **node2,
    struct spamm_ll_t *list);

#endif
