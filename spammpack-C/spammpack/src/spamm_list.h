/** @file */

#ifndef __SPAMM_LIST_H
#define __SPAMM_LIST_H

struct spamm_list_t;

struct spamm_list_t *
spamm_list_new (const unsigned int length);

void
spamm_list_delete (struct spamm_list_t **list);

int
spamm_list_compare_int (const unsigned int a, const unsigned int b, void *user_data);

void
spamm_list_sort (struct spamm_list_t *list,
    int (*compare) (const unsigned int, const unsigned int, void *),
    void *user_data);

unsigned int
spamm_list_length (struct spamm_list_t *list);

unsigned int
spamm_list_get (struct spamm_list_t *list, const unsigned int i);

void
spamm_list_set (struct spamm_list_t *list, const unsigned int i, const unsigned int Ai);

#endif
