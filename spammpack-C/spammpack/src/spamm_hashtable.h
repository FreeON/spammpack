/** @file */

#ifndef __SPAMM_HASHTABLE_H
#define __SPAMM_HASHTABLE_H

#include "spamm_list.h"

#include <stdint.h>

struct spamm_hashtable_t;

uint32_t
spamm_hashtable_hash (const uint32_t key);

uint32_t
spamm_hashtable_hash_Jenkins (const uint32_t key);

uint32_t
spamm_hashtable_hash_MurmurHash3 (const uint32_t key);

uint32_t
spamm_hashtable_hash_direct (const uint32_t key);

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
spamm_hashtable_keys (const struct spamm_hashtable_t *hashtable);

unsigned int
spamm_hashtable_get_number_buckets (const struct spamm_hashtable_t *hashtable);

unsigned int
spamm_hashtable_get_total_distance (const struct spamm_hashtable_t *hashtable);

unsigned int
spamm_hashtable_get_number_keys (const struct spamm_hashtable_t *hashtable);

unsigned int
spamm_hashtable_memory (const struct spamm_hashtable_t *hashtable);

#endif
