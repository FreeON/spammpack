#include <spamm_hashtable.h>
#include <stdio.h>
#include <stdlib.h>

#define NUMBER_KEYS 5000

int
main ()
{
  int result = 0;
  int i;
  int *value;
  struct spamm_hashtable_t *hashtable;

  hashtable = spamm_hashtable_new();

  /* Insert keys. */
  for(i = 0; i < NUMBER_KEYS; i++)
  {
    value = malloc(sizeof(int));
    *value = i;
    spamm_hashtable_insert(hashtable, i, value);
  }

  /* Test whether all those stored keys are there. */
  for(i = 0; i < NUMBER_KEYS; i++)
  {
    value = spamm_hashtable_lookup(hashtable, i);
    if(value == NULL)
    {
      printf("found NULL value for key %u\n", i);
      return 1;
    }

    else if(*((int*) value) != i)
    {
      printf("incorrect value for key %u, value = %u\n", i, *value);
      return 2;
    }

    if((value = spamm_hashtable_remove(hashtable, i)) == NULL)
    {
      printf("could not remove key %u\n", i);
      return 3;
    }

    else { free(value); }
  }

  spamm_hashtable_delete(&hashtable);

  return result;
}
