#include "spamm.h"
#include "spamm_hashtable.h"
#include "spamm_list.h"
#include "spamm_types.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

/** An empty key. */
#define SPAMM_KEY_EMPTY 0

/** A deleted key. */
#define SPAMM_KEY_DELETED 1

/* Define false. */
#define SPAMM_FALSE 0

/** Define true. */
#define SPAMM_TRUE 1

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

#ifdef SPAMM_USE_HASHTABLE
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
#else
  /** The stored values. */
  void **value;
#endif
};

/** A hash function. This is a wrapper for the actual hash function used.
 *
 * @param key The key to hash.
 *
 * @return The hashed key.
 */
inline uint32_t
spamm_hashtable_hash (const uint32_t key)
{
  return spamm_hashtable_hash_direct(key);
  //return spamm_hashtable_hash_Jenkins(key);
  //return spamm_hashtable_hash_MurmurHash3(key);
}

/** A hash function. This is taken the hashword() function from Bob Jenkins'
 * webpage,
 *
 * http://burtleburtle.net/bob/c/lookup3.c
 *
 * @param key The key to hash.
 *
 * @return The hashed key.
 */
uint32_t
spamm_hashtable_hash_Jenkins (const uint32_t key)
{
  static const uint32_t seed = 0xC80748E8;
  unsigned int a, b, c;

  /* Set up the internal state */
  a = b = c = 0xdeadbeef + (((uint32_t) 4) + seed);

  /* Rotate. */
  a += key;
  c ^= b; c -= ((b << 14) | (b >> (32 - 14)));
  a ^= c; a -= ((c << 11) | (c >> (32 - 11)));
  b ^= a; b -= ((a << 25) | (a >> (32 - 25)));
  c ^= b; c -= ((b << 16) | (b >> (32 - 16)));
  a ^= c; a -= ((c <<  4) | (c >> (32 -  4)));
  b ^= a; b -= ((a << 14) | (a >> (32 - 14)));
  c ^= b; c -= ((b << 24) | (b >> (32 - 24)));

  /* Return the result. */
  return c;
}

/** A really simple hash function. We assume here that we store a full matrix,
 * i.e. all linear indices. This means that we can pass the key right through
 * and don't have to do anything.
 *
 * @param key The key to hash.
 *
 * @return The hashed key.
 */
inline uint32_t
spamm_hashtable_hash_direct (const uint32_t key)
{
  return key;
}

/** A hash function. This one is taken from the MurmurHash3_x86_32() function
 * of MurmuHash3.cpp from Austin Appleby.
 *
 * http://code.google.com/p/smhasher/
 *
 * @param key The key (32-bit).
 * @param seed The seed.
 *
 * @return The hashed key.
 */
uint32_t
spamm_hashtable_hash_MurmurHash3 (const uint32_t key)
{
#define HASH_MURMURHASH3_SEED 1

  const uint8_t *data = (const uint8_t*) &key;
  uint32_t h1 = HASH_MURMURHASH3_SEED;
  uint32_t c1 = 0xcc9e2d51;
  uint32_t c2 = 0x1b873593;
  uint32_t k1 = 0;

  k1 ^= data[0];
  k1 *= c1;
  k1  = (k1 << 15) | (k1 >> (32 - 15));;
  k1 *= c2;
  h1 ^= k1;

  h1 ^= 1;

  h1 ^= h1 >> 16;
  h1 *= 0x85ebca6b;
  h1 ^= h1 >> 13;
  h1 *= 0xc2b2ae35;
  h1 ^= h1 >> 16;

  return h1;
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
  struct spamm_hashtable_t *hashtable;
#ifdef SPAMM_USE_HASHTABLE
  short i;

  /* We should have at least 8 buckets. */
  if (number_buckets < 8) { number_buckets = 8; }

  /* The number of buckets should be a power of 2. */
  number_buckets--;
  for (i = 1; i < sizeof(unsigned int)*8-1; i++)
  {
    number_buckets |= number_buckets >> i;
  }
  number_buckets++;

  hashtable = calloc(1, sizeof(struct spamm_hashtable_t));
  hashtable->number_buckets = number_buckets;
  hashtable->data = calloc(hashtable->number_buckets, sizeof(struct spamm_hashtable_bucket_t));
#else
  hashtable = calloc(1, sizeof(struct spamm_hashtable_t));
  hashtable->number_buckets = number_buckets;
  hashtable->value = calloc(hashtable->number_buckets, sizeof(void*));
#endif

  return hashtable;
}

/** Create a new hashtable.
 *
 * @return The newly allocated hashtable.
 */
struct spamm_hashtable_t *
spamm_hashtable_new ()
{
  return spamm_hashtable_new_sized(4);
}

/** Destroy a hashtable. Note that the values stored are not free'ed.
 *
 * @param hashtable The hashtable to free.
 */
void
spamm_hashtable_delete (struct spamm_hashtable_t **hashtable)
{
#ifdef SPAMM_USE_HASHTABLE
  free((*hashtable)->data);
#else
  free((*hashtable)->value);
#endif

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
#ifdef SPAMM_USE_HASHTABLE
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
    if (old_hashtable->has_special_key[i] == SPAMM_TRUE)
    {
      new_hashtable->has_special_key[i] = SPAMM_TRUE;
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
#else
  unsigned int i;
  void **value = calloc(number_buckets, sizeof(void*));

  /* Copy old data. */
  for (i = 0; i < hashtable->number_buckets; i++)
  {
    value[i] = hashtable->value[i];
  }

  /* Free old data array. */
  free(hashtable->value);

  /* Copy data into hashtable. */
  hashtable->number_buckets = number_buckets;
  hashtable->value = value;
#endif
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
#ifdef SPAMM_USE_HASHTABLE
  unsigned int keyhash = spamm_hashtable_hash(key);
  unsigned int bucket_index = keyhash & (hashtable->number_buckets-1);

  if (key == SPAMM_KEY_EMPTY || key == SPAMM_KEY_DELETED)
  {
    /* The special keys SPAMM_KEY_EMPTY and SPAMM_KEY_DELETED have to be
     * treated first because the usual key hashing will not work properly on
     * them.
     */
    hashtable->has_special_key[key] = SPAMM_TRUE;
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
    spamm_hashtable_insert(hashtable, key, value);
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
#else
  if (key < hashtable->number_buckets)
  {
    hashtable->value[key] = value;
  }

  else
  {
    spamm_hashtable_rehash(hashtable, hashtable->number_buckets << 1);
    spamm_hashtable_insert(hashtable, key, value);
  }
#endif
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
#ifdef SPAMM_USE_HASHTABLE
  unsigned int keyhash = spamm_hashtable_hash(key);
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
#endif

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
#ifdef SPAMM_USE_HASHTABLE
  void *value = NULL;
  struct spamm_hashtable_bucket_t *bucket;

  if (key == SPAMM_KEY_EMPTY || key == SPAMM_KEY_DELETED)
  {
    /* The special keys SPAMM_KEY_EMPTY and SPAMM_KEY_DELETED have to be
     * treated first because the usual key hashing will not work properly on
     * them.
     */
    if (hashtable->has_special_key[key] == SPAMM_TRUE)
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
#else
  if (key < hashtable->number_buckets)
  {
    return hashtable->value[key];
  }

  else
  {
    spamm_hashtable_rehash(hashtable, hashtable->number_buckets << 1);
    return spamm_hashtable_lookup(hashtable, key);
  }
#endif
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
#ifdef SPAMM_USE_HASHTABLE
  struct spamm_hashtable_bucket_t *bucket = NULL;

  if (key == SPAMM_KEY_EMPTY || key == SPAMM_KEY_DELETED)
  {
    if (hashtable->has_special_key[key] == SPAMM_TRUE)
    {
      hashtable->has_special_key[key] = SPAMM_FALSE;
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
#else
  if (key < hashtable->number_buckets)
  {
    result = hashtable->value[key];
    hashtable->value[key] = NULL;
  }

  else
  {
    printf("[spamm_hashtable_remove] key out of bounds\n");
    exit(1);
  }
#endif

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
    void (*func) (unsigned int, void *, void *), void *user_data)
{
  unsigned int i;

#ifdef SPAMM_USE_HASHTABLE
  /* Loop over all buckets and call func() for each key. */
  for (i = 0; i < 2; i++)
  {
    if (hashtable->has_special_key[i] == SPAMM_TRUE)
    {
      func(i, hashtable->special_value[i], user_data);
    }
  }

  for (i = 0; i < hashtable->number_buckets; i++)
  {
    if (hashtable->data[i].key != SPAMM_KEY_EMPTY &&
        hashtable->data[i].key != SPAMM_KEY_DELETED)
    {
      func(hashtable->data[i].key, hashtable->data[i].value, user_data);
    }
  }
#else
  for (i = 0; i < hashtable->number_buckets; i++)
  {
    if (hashtable->value[i] != NULL)
    {
      func(i, hashtable->value[i], user_data);
    }
  }
#endif
}

/** Get a list of matrix indices and a list of norms.
 *
 * @param hashtable The hashtable.
 *
 * @return list An empty list for the matrix indices. The list will be
 * allocated in this function and needs to be free()'ed by the caller when not
 * needed anymore.
 */
struct spamm_list_t *
spamm_hashtable_keys (const struct spamm_hashtable_t *hashtable)
{
  struct spamm_list_t *list;
  unsigned int i;
  unsigned int list_index;

#ifdef SPAMM_USE_HASHTABLE
  /* Allocate new lists. */
  list = spamm_list_new((hashtable->has_special_key[0] == SPAMM_TRUE ? 1 : 0)
      + (hashtable->has_special_key[1] == SPAMM_TRUE ? 1 : 0)
      + hashtable->number_stored_keys);

  /* Populate lists. */
  list_index = 0;

  for (i = 0; i < 2; i++)
  {
    if (hashtable->has_special_key[i] == SPAMM_TRUE)
    {
      spamm_list_set(list, list_index, i, ((struct spamm_data_t*) hashtable->special_value[i])->node_norm);
      list_index++;
    }
  }

  for (i = 0; i < hashtable->number_buckets; i++)
  {
    if (hashtable->data[i].key != SPAMM_KEY_EMPTY &&
        hashtable->data[i].key != SPAMM_KEY_DELETED)
    {
      spamm_list_set(list, list_index, hashtable->data[i].key, ((struct spamm_data_t*) hashtable->data[i].value)->node_norm);
      list_index++;
    }
  }

  if (list_index != spamm_list_length(list))
  {
    printf("error copying indices\n");
    exit(1);
  }
#else
  unsigned int number_keys;

  /* Count number of keys (the entries that are not NULL).
   */
  number_keys = 0;
  for (i = 0; i < hashtable->number_buckets; i++)
  {
    if (hashtable->value[i] != NULL)
    {
      number_keys++;
    }
  }

  /* Copy the existing keys into a list. */
  list = spamm_list_new(number_keys);
  list_index = 0;
  for (i = 0; i < hashtable->number_buckets; i++)
  {
    if (hashtable->value[i] != NULL)
    {
      spamm_list_set(list, list_index, i, ((struct spamm_data_t*) hashtable->value[i])->node_norm);
      list_index++;
    }
  }
#endif

  return list;
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
#ifdef SPAMM_USE_HASHTABLE
  return hashtable->total_distance;
#else
  return 0;
#endif
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
#ifdef SPAMM_USE_HASHTABLE
  short i;
  unsigned int result = hashtable->number_stored_keys;

  for (i = 0; i < 2; i++)
  {
    if (hashtable->has_special_key[i] == SPAMM_TRUE)
    {
      result++;
    }
  }

  return hashtable->number_stored_keys;
#else
  return hashtable->number_buckets;
#endif
}

/** Return the memory used in a hashtable. This does <em>not</em> include the
 * memory consumed by the values.
 *
 * @param hashtable The hashtable.
 *
 * @return The number of bytes consumed by the hashtable.
 */
unsigned int
spamm_hashtable_memory (const struct spamm_hashtable_t *hashtable)
{
  unsigned int total = 0;

  total = sizeof(struct spamm_hashtable_t);
  total += hashtable->number_buckets*sizeof(struct spamm_hashtable_bucket_t);

  return total;
}
