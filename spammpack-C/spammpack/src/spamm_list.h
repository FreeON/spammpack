/** @file */

#ifndef __SPAMM_LIST_H
#define __SPAMM_LIST_H

struct spamm_list_t;

/** The comparison function of 2 indices for the sort. */
typedef int (*spamm_list_compare_index_function) (const unsigned int, const unsigned int);

struct spamm_list_t *
spamm_list_new (const unsigned int length);

void
spamm_list_delete (struct spamm_list_t **list);

int
spamm_list_compare_int (const unsigned int a, const unsigned int b);

int
spamm_list_compare_index_row (const unsigned int a, const unsigned int b);

int
spamm_list_compare_index_column (const unsigned int a, const unsigned int b);

void
spamm_list_sort_index (struct spamm_list_t *list, spamm_list_compare_index_function compare);

void
spamm_list_sort_norm (struct spamm_list_t *list, const unsigned int left, const unsigned int right);

unsigned int
spamm_list_length (struct spamm_list_t *list);

unsigned int
spamm_list_get_index (struct spamm_list_t *list, const unsigned int i);

float
spamm_list_get_norm (struct spamm_list_t *list, const unsigned int i);

void
spamm_list_set (struct spamm_list_t *list, const unsigned int i, const unsigned int index, const float norm);

void
spamm_list_print (const struct spamm_list_t *list);

#endif
