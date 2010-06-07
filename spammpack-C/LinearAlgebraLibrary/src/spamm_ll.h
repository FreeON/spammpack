/** @file */

#if defined(__SPAMM_LL_H)
#warn Already included spamm_ll.h
#else

/** Define in case spamm_ll.h has been included. */
#define __SPAMM_LL_H 1

#include "config.h"

/** \page page_ll The SpAMM linked list.
 *
 * Some functions for dealing with linked lists. Details can be found in
 * spamm_ll.h.
 *
 * \author Nicolas Bock <nicolasbock@gmail.com>
 * \author Matt Challacombe <matt.challacombe@gmail.com>.
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

void
spamm_ll_initialize (struct spamm_ll_t *list);

void
spamm_ll_initialize_node (struct spamm_ll_node_t **node);

void
spamm_ll_delete (struct spamm_ll_t *list);

void
spamm_ll_append (void *data, struct spamm_ll_t *list);

struct spamm_ll_node_t *
spamm_ll_get (const unsigned int i, const struct spamm_ll_t *list);

void
spamm_ll_swap (struct spamm_ll_node_t **node1, struct spamm_ll_node_t **node2,
    struct spamm_ll_t *list);

void
spamm_ll_sort (int (*compare) (const void *data1, const void *data2), struct spamm_ll_t *list);

void
spamm_ll_print_node (const struct spamm_ll_node_t *node);

void
spamm_ll_print (const struct spamm_ll_t *list);

#endif
