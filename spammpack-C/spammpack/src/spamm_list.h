/** @file */

#ifndef __SPAMM_LIST_H
#define __SPAMM_LIST_H

struct spamm_list_t;

typedef int (*spamm_list_compare_function) (const unsigned int, const unsigned int, const float, const float);

struct spamm_list_t *
spamm_list_new (const unsigned int length);

void
spamm_list_delete (struct spamm_list_t **list);

int
spamm_list_compare_int (const unsigned int a, const unsigned int b, const float a_norm, const float b_norm);

void
spamm_list_sort (
    struct spamm_list_t *list,
    float *node_norm,
    spamm_list_compare_function compare);

unsigned int
spamm_list_length (struct spamm_list_t *list);

unsigned int
spamm_list_get (struct spamm_list_t *list, const unsigned int i);

unsigned int*
spamm_list_get_data (struct spamm_list_t *list);

void
spamm_list_set (struct spamm_list_t *list, const unsigned int i, const unsigned int Ai);

#endif
