/** @file */

#ifndef __SPAMM_HASHTABLE_H
#define __SPAMM_HASHTABLE_H

struct spamm_hashtable_t *
spamm_hashtable_new ();

void
spamm_hashtable_delete (struct spamm_hashtable_t **hashtable);

void
spamm_hashtable_insert (struct spamm_hashtable_t *hashtable,
    unsigned int key, void *value);

void *
spamm_hashtable_lookup (const struct spamm_hashtable_t *hashtable,
    const unsigned int key);

#endif
