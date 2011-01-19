#include "spamm.h"
#include <stdio.h>
#include <stdlib.h>

#define N 10000

int
main ()
{
  int result = 0;
  unsigned int i;
  unsigned int last_i;
  unsigned int value, last_value;
  struct spamm_list_t *list;

  list = spamm_list_new(N);
  for (i = 0; i < N; i++)
  {
    spamm_list_set(list, i, rand());
  }

  spamm_list_sort(list, spamm_list_compare_int, NULL);

  last_i = 0;
  last_value = spamm_list_get(list, 0);
  for (i = 0; i < N; i++)
  {
    value = spamm_list_get(list, i);
    if (value > last_value)
    {
      last_value = value;
      last_i = i;
    }

    else if (value < last_value)
    {
      result = -1;
      printf("failed to sort, list[%u] = %i, list[%u] = %i...\n", i, spamm_list_get(list, i), last_i, last_value);
      break;
    }
  }

  spamm_list_delete(&list);

  return result;
}
