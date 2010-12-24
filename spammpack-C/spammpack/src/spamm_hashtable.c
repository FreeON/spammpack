#include "spamm_hashtable.h"
#include "spamm_list.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/** An empty key. */
#define SPAMM_KEY_EMPTY 0

/** A deleted key. */
#define SPAMM_KEY_DELETED 1

/** The grow threshold. */
#define SPAMM_GROW_THRESHOLD 0.65

/** The shrink threshold. */
#define SPAMM_SHRINK_THRESHOLD 0.30

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

  /** A measure of the number of collisions in the hashtable. We measure the
   * "distance" of the inserted key to where it originally was going to be
   * inserted. Since we deal with collision by linear chaining, a distance
   * greater than zero indicates a deviation of the scaling from O(1) to O(N).
   */
  unsigned int total_distance;

  /** The buckets. */
  struct spamm_hashtable_bucket_t *data;
};

/** A utility function for Bob Jenkins' hash function. */
#define rot(x,k) (((x) << (k)) | ((x) >> (32-(k))))

/** A hash function. This is taken the hashword() function from Bob Jenkins'
 * webpage,
 *
 * http://burtleburtle.net/bob/c/lookup3.c
 *
 * @param key The key to hash.
 *
 * @return The hashed key.
 */
unsigned int
spamm_hashtable_hash_1 (const unsigned int key)
{
  unsigned int a, b, c;

  /* Set up the internal state */
  a = b = c = 0xdeadbeef + (((uint32_t) 1) << 2);

  /* Rotate. */
  a += key;
  c ^= b; c -= rot(b,14);
  a ^= c; a -= rot(c,11);
  b ^= a; b -= rot(a,25);
  c ^= b; c -= rot(b,16);
  a ^= c; a -= rot(c,4);
  b ^= a; b -= rot(a,14);
  c ^= b; c -= rot(b,24);

  /* Return the result. */
  return c;
}

/** Create a new hashtable.
 *
 * @param number_buckets The number of buckets to create.
 *
 * @return The newly allocated hashtable.
 */
struct spamm_hashtable_t *
spamm_hashtable_new_sized (unsigned int number_buckets)
{
  short i;
  struct spamm_hashtable_t *hashtable;

  /* We should have at least 8 buckets. */
  if (number_buckets < 8) { number_buckets = 8; }

  /* The number of buckets should be a power of 2. */
  number_buckets--;
  for (i = 1; i < sizeof(unsigned int)*8-1; i++)
  {
    number_buckets |= number_buckets >> i;
  }
  number_buckets++;

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
  return spamm_hashtable_new_sized(8);
}

/** Destroy a hashtable. Note that the values stored are freed with free().
 *
 * @param hashtable The hashtable to free.
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

  new_hashtable = spamm_hashtable_new_sized(number_buckets);

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
  old_hashtable->total_distance = new_hashtable->total_distance;
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
  unsigned int keyhash = spamm_hashtable_hash_1(key);
  unsigned int bucket_index = keyhash & (hashtable->number_buckets-1);

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
      hashtable->total_distance++;
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
    while (hashtable->data[bucket_index].key != SPAMM_KEY_EMPTY);
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

/** Lookup a key in a hashtable and return a bucket item.
 *
 * @param hashtable The hashtable.
 * @param key The key to lookup.
 *
 * @return NULL in case the key is not in the table, or if the value stored
 * there is NULL. Any value other than NULL indicates that they key was found,
 * and the pointer returned is the pointer stored with this key.
 */
struct spamm_hashtable_bucket_t *
spamm_hashtable_lookup_bucket (struct spamm_hashtable_t *hashtable,
    const unsigned int key)
{
  unsigned int keyhash = spamm_hashtable_hash_1(key);
  unsigned int bucket_index = keyhash & (hashtable->number_buckets-1);

  if (hashtable->data[bucket_index].key == key)
  {
    return &hashtable->data[bucket_index];
  }

  else if (hashtable->data[bucket_index].key == SPAMM_KEY_EMPTY)
  {
    return NULL;
  }

  else
  {
    /* Linear chaining, i.e. increment bucket index until we either find the
     * key or an empty bucket. */
    do
    {
      bucket_index = (bucket_index+1) & (hashtable->number_buckets-1);
      if (hashtable->data[bucket_index].key == key)
      {
        return &hashtable->data[bucket_index];
      }
    }
    while (hashtable->data[bucket_index].key != SPAMM_KEY_EMPTY);
  }

  return NULL;
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
  struct spamm_hashtable_bucket_t *bucket;

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

  else
  {
    bucket = spamm_hashtable_lookup_bucket(hashtable, key);
    if (bucket != NULL)
    {
      return bucket->value;
    }

    else
    {
      return NULL;
    }
  }
}

/** Remove a key from a hashtable.
 *
 * @param hashtable The hashtable.
 * @param key The key to remove.
 *
 * @return A pointer to the removed value. If this pointer is NULL then either
 * they key was not found or the value stored was actually NULL.
 */
void *
spamm_hashtable_remove (struct spamm_hashtable_t *hashtable,
    const unsigned int key)
{
  void *result = NULL;
  struct spamm_hashtable_bucket_t *bucket = NULL;

  if (key == SPAMM_KEY_EMPTY || key == SPAMM_KEY_DELETED)
  {
    if (hashtable->has_special_key[key] == 1)
    {
      hashtable->has_special_key[key] = 0;
      result = hashtable->special_value[key];
      hashtable->special_value[key] = NULL;
    }
  }

  else
  {
    bucket = spamm_hashtable_lookup_bucket(hashtable, key);
    if (bucket != NULL)
    {
      bucket->key = SPAMM_KEY_DELETED;
      result = bucket->value;
      bucket->value = NULL;

      hashtable->number_stored_keys--;
      hashtable->number_deleted_keys++;

      if (hashtable->number_stored_keys < SPAMM_SHRINK_THRESHOLD*hashtable->number_buckets)
      {
        spamm_hashtable_rehash(hashtable, hashtable->number_buckets >> 1);
      }
    }
  }

  return result;
}

/** Iterate over keys and call a function for each key.
 *
 * @param hashtable The hashtable.
 * @param func The function to execute for each of the keys.
 * @param user_data A pointer that will be passed to func.
 */
void
spamm_hashtable_foreach (struct spamm_hashtable_t *hashtable,
    void (*func) (void *, void *, void *), void *user_data)
{
}

/** Copy all keys into array. The returned pointer points to an array that is
 * allocated in this function with malloc(). The caller needs to take care of
 * calling free() on it when the array is not needed anymore.
 *
 * @param hashtable The hashtable.
 *
 * @return A pointer to an array that contains all of the keys.
 */
struct spamm_list_t *
spamm_hashtable_keys (struct spamm_hashtable_t *hashtable)
{
}

/** Get the number of buckets in this hashtable.
 *
 * @param hashtable The hashtable.
 *
 * @return The number of buckets.
 */
unsigned int
spamm_hashtable_get_number_buckets (const struct spamm_hashtable_t *hashtable)
{
  return hashtable->number_buckets;
}

/** Get the total distance of the keys in the hashtable.
 *
 * @param hashtable The hashtable.
 *
 * @return The total distance.
 */
unsigned int
spamm_hashtable_get_total_distance (const struct spamm_hashtable_t *hashtable)
{
  return hashtable->total_distance;
}

/** Get the number of stored keys.
 *
 * @param hashtable The hashtable.
 *
 * @return The number of stored keys.
 */
unsigned int
spamm_hashtable_get_number_keys (const struct spamm_hashtable_t *hashtable)
{
  short i;
  unsigned int result = hashtable->number_stored_keys;

  for (i = 0; i < 2; i++)
  {
    if (hashtable->has_special_key[i] == 1)
    {
      result++;
    }
  }

  return hashtable->number_stored_keys;
}
