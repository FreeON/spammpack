/** @file */

#include <glib.h>

#ifndef __SPAMM_HASHTABLE_H
#define __SPAMM_HASHTABLE_H

gboolean
spamm_hash_uint_equal (gconstpointer a, gconstpointer b);

unsigned int
spamm_hashtable_hash_1 (const unsigned int key);

guint
spamm_hashtable_hash_1_glib (gconstpointer key);

struct spamm_hashtable_t *
spamm_hashtable_new ();

struct spamm_hashtable_t *
spamm_hashtable_new_sized (const unsigned int number_buckets);

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

unsigned int
spamm_hashtable_get_number_buckets (const struct spamm_hashtable_t *hashtable);

unsigned int
spamm_hashtable_get_total_distance (const struct spamm_hashtable_t *hashtable);

unsigned int
spamm_hashtable_get_number_keys (const struct spamm_hashtable_t *hashtable);

#endif
