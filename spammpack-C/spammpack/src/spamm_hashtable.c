#include "spamm.h"
#include <stdio.h>
#include <stdlib.h>

/** Compare to unsigned int values. This is a helper function for the
 * hashtables in dealing with the linear indices.
 *
 * @param a A pointer to the first index.
 * @param b A pointer to the second index.
 *
 * @return TRUE if the two values are equal, FALSE if not.
 */
gboolean
spamm_hash_uint_equal (gconstpointer a, gconstpointer b)
{
  const unsigned int *value_1 = a;
  const unsigned int *value_2 = b;

  if ((*value_1) == (*value_2)) { return TRUE; }
  else { return FALSE; }
}

/** An empty key. */
#define SPAMM_KEY_EMPTY 0

/** A deleted key. */
#define SPAMM_KEY_DELETED 1

/** The grow threshold. */
#define SPAMM_GROW_THRESHOLD 0.65

/** A hashtable bucket.
 */
struct spamm_hashtable_bucket_t
{
  /** The key stored in this bucket. */
  unsigned int key;

  /** The value stored in this bucket. */
  void *value;
};

/** A hashtable.
 */
struct spamm_hashtable_t
{
  /** The number of buckets in this hashtable. */
  unsigned int number_buckets;

  /** Flag for whether the special keys SPAMM_KEY_EMPTY and SPAMM_KEY_DELETE
   * are stored in this hashtable.
   */
  short has_special_key[2];

  /** Values stored for the special keys. */
  void *special_value[2];

  /** Number of deleted keys. */
  unsigned int number_deleted_keys;

  /** Number of non-empty keys. */
  unsigned int number_stored_keys;

  /** The buckets. */
  struct spamm_hashtable_bucket_t *data;
};

/** Create a new hashtable (internal version).
 *
 * @param number_buckets The number of buckets to create.
 *
 * @return The newly allocated hashtable.
 */
struct spamm_hashtable_t *
spamm_hashtable_new_internal (const unsigned int number_buckets)
{
  struct spamm_hashtable_t *hashtable;

  hashtable = (struct spamm_hashtable_t*) calloc(sizeof(struct spamm_hashtable_t), 1);
  hashtable->number_buckets = number_buckets;
  hashtable->data = calloc(sizeof(struct spamm_hashtable_bucket_t), hashtable->number_buckets);

  return hashtable;
}

/** Create a new hashtable.
 *
 * @return The newly allocated hashtable.
 */
struct spamm_hashtable_t *
spamm_hashtable_new ()
{
  return spamm_hashtable_new_internal(8);
}

/** Destroy a hashtable. Note that the values stored are freed with free().
 *
 * @param hastable The hashtable to free.
 */
void
spamm_hashtable_delete (struct spamm_hashtable_t **hashtable)
{
  unsigned int i;

  for (i = 0; i < (*hashtable)->number_buckets; i++)
  {
    if ((*hashtable)->data[i].key != SPAMM_KEY_EMPTY &&
        (*hashtable)->data[i].key != SPAMM_KEY_DELETED)
    {
      free((*hashtable)->data[i].value);
    }
  }

  for (i = 0; i < 2; i++)
  {
    if ((*hashtable)->has_special_key[i] == 1)
    {
      free((*hashtable)->special_value[i]);
    }
  }

  free((*hashtable)->data);
  free(*hashtable);
  *hashtable = NULL;
}

void
spamm_hashtable_insert (struct spamm_hashtable_t *hashtable,
    unsigned int key, void *value);

/** Rehash a hashtable to resize the number of buckets.
 *
 * @param hashtable The hashtable.
 * @param number_buckets The number of buckets to resize to.
 */
void
spamm_hashtable_rehash (struct spamm_hashtable_t *hashtable,
    const unsigned int number_buckets)
{
  unsigned int i;
  struct spamm_hashtable_t *old_hashtable = hashtable;
  struct spamm_hashtable_t *new_hashtable;

  new_hashtable = spamm_hashtable_new_internal(number_buckets);

  /* Re-insert all keys into new hashtable. */
  for (i = 0; i < old_hashtable->number_buckets; i++)
  {
    if (old_hashtable->data[i].key != SPAMM_KEY_EMPTY &&
        old_hashtable->data[i].key != SPAMM_KEY_DELETED)
    {
      spamm_hashtable_insert(new_hashtable, old_hashtable->data[i].key, old_hashtable->data[i].value);
    }
  }

  /* Sanity check. */
  if (old_hashtable->number_stored_keys != new_hashtable->number_stored_keys)
  {
    printf("error in hashtable rehash... had %u keys in old hashtable, have %u keys in new hashtable\n",
        old_hashtable->number_stored_keys, new_hashtable->number_stored_keys);
    exit(1);
  }

  /* Copy the 2 special keys. */
  for (i = 0; i < 2; i++)
  {
    if (old_hashtable->has_special_key[i] == 1)
    {
      new_hashtable->has_special_key[i] = 1;
      new_hashtable->special_value[i] = old_hashtable->special_value[i];
    }
  }

  /* Delete buckets of old hashtable. */
  free(old_hashtable->data);

  /* Copy all information of new hashtable into old one. */
  old_hashtable->number_buckets = new_hashtable->number_buckets;
  old_hashtable->number_deleted_keys = new_hashtable->number_deleted_keys;
  old_hashtable->number_stored_keys = new_hashtable->number_stored_keys;
  old_hashtable->data = new_hashtable->data;

  /* Free the new hashtable. */
  free(new_hashtable);
}

/** Insert a key/value pair into the hashtable. If the key to be inserted
 * already exists in the hashtable, the value stored is overwritten with the
 * new value.
 *
 * @param hashtable The hashtable.
 * @param key The key to insert.
 * @param value The value to insert.
 */
void
spamm_hashtable_insert (struct spamm_hashtable_t *hashtable,
    unsigned int key, void *value)
{
  unsigned int bucket_index = key & (hashtable->number_buckets-1);

  if (key == SPAMM_KEY_EMPTY || key == SPAMM_KEY_DELETED)
  {
    /* The special keys SPAMM_KEY_EMPTY and SPAMM_KEY_DELETED have to be
     * treated first because the usual key hashing will not work properly on
     * them.
     */
    hashtable->has_special_key[key] = 1;
    hashtable->special_value[key] = value;
    return;
  }

  else if (hashtable->data[bucket_index].key == key)
  {
    /* Store value in this bucket. This overwrites any existing value. */
    hashtable->data[bucket_index].value = value;
    return;
  }

  else if (hashtable->data[bucket_index].key != SPAMM_KEY_EMPTY &&
      hashtable->data[bucket_index].key != SPAMM_KEY_DELETED)
  {
    /* Linear chaining, i.e. increment bucket index until we either find the
     * key or an empty bucket. */
    do
    {
      bucket_index = (bucket_index+1) & (hashtable->number_buckets-1);
      if (hashtable->data[bucket_index].key == key)
      {
        hashtable->data[bucket_index].value = value;
        return;
      }

      else if (hashtable->data[bucket_index].key == SPAMM_KEY_DELETED)
      {
        break;
      }
    }
    while (hashtable->data[bucket_index].key > SPAMM_KEY_EMPTY);
  }

  if (hashtable->data[bucket_index].key == SPAMM_KEY_EMPTY &&
      hashtable->number_stored_keys >= SPAMM_GROW_THRESHOLD*hashtable->number_buckets)
  {
    spamm_hashtable_rehash(hashtable, hashtable->number_buckets << 1);
    return spamm_hashtable_insert(hashtable, key, value);
  }

  else
  {
    if (hashtable->data[bucket_index].key == SPAMM_KEY_DELETED)
    {
      hashtable->number_deleted_keys--;
    }

    hashtable->data[bucket_index].key = key;
    hashtable->data[bucket_index].value = value;
    hashtable->number_stored_keys++;
    return;
  }
}

/** Lookup a key in the hashtable.
 *
 * @param hashtable The hashtable.
 * @param key The key to lookup.
 *
 * @return NULL in case the key is not in the table, or if the value stored
 * there is NULL. Any value other than NULL indicates that they key was found,
 * and the pointer returned is the pointer stored with this key.
 */
void *
spamm_hashtable_lookup (struct spamm_hashtable_t *hashtable,
    const unsigned int key)
{
  void *value = NULL;
  unsigned int bucket_index = key & (hashtable->number_buckets-1);

  if (key == SPAMM_KEY_EMPTY || key == SPAMM_KEY_DELETED)
  {
    /* The special keys SPAMM_KEY_EMPTY and SPAMM_KEY_DELETED have to be
     * treated first because the usual key hashing will not work properly on
     * them.
     */
    if (hashtable->has_special_key[key] == 1)
    {
      value = hashtable->special_value[key];
    }
    return value;
  }

  else if (hashtable->data[bucket_index].key == key)
  {
    /* Store value in this bucket. This overwrites any existing value. */
    return hashtable->data[bucket_index].value;
  }

  else if (hashtable->data[bucket_index].key != SPAMM_KEY_EMPTY &&
      hashtable->data[bucket_index].key != SPAMM_KEY_DELETED)
  {
    /* Linear chaining, i.e. increment bucket index until we either find the
     * key or an empty bucket. */
    do
    {
      bucket_index = (bucket_index+1) & (hashtable->number_buckets-1);
      if (hashtable->data[bucket_index].key == key)
      {
        return hashtable->data[bucket_index].value;
      }

      else if (hashtable->data[bucket_index].key == SPAMM_KEY_DELETED)
      {
        break;
      }
    }
    while (hashtable->data[bucket_index].key > SPAMM_KEY_EMPTY);
  }
}
