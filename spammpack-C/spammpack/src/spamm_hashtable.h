/** @file */

#ifndef __SPAMM_HASHTABLE_H
#define __SPAMM_HASHTABLE_H

#include "spamm_list.h"

struct spamm_hashtable_t;

unsigned int
spamm_hashtable_hash_1 (const unsigned int key);

struct spamm_hashtable_t *
spamm_hashtable_new ();

struct spamm_hashtable_t *
spamm_hashtable_new_sized (const unsigned int number_buckets);

void
spamm_hashtable_foreach (struct spamm_hashtable_t *hashtable,
    void (*func) (unsigned int, void *, void *), void *user_data);

void
spamm_hashtable_delete (struct spamm_hashtable_t **hashtable);

void
spamm_hashtable_insert (struct spamm_hashtable_t *hashtable,
    unsigned int key, void *value);

void *
spamm_hashtable_lookup (struct spamm_hashtable_t *hashtable,
    const unsigned int key);

void *
spamm_hashtable_remove (struct spamm_hashtable_t *hashtable,
    const unsigned int key);

struct spamm_list_t *
spamm_hashtable_keys (struct spamm_hashtable_t *hashtable);

unsigned int
spamm_hashtable_get_number_buckets (const struct spamm_hashtable_t *hashtable);

unsigned int
spamm_hashtable_get_total_distance (const struct spamm_hashtable_t *hashtable);

unsigned int
spamm_hashtable_get_number_keys (const struct spamm_hashtable_t *hashtable);

unsigned int
spamm_hashtable_memory (const struct spamm_hashtable_t *hashtable);

#endif
