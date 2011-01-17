#include "spamm.h"
#include <stdio.h>
#include <stdlib.h>

#define N 100

int
main ()
{
  int result = 0;
  int i;
  struct spamm_list_t *list;

  list = spamm_list_new(N);
  for (i = 0; i < N; i++)
  {
    spamm_list_set(list, i, N-i-1);
  }

  spamm_list_sort(list, spamm_list_compare_int, NULL);

  for (i = 0; i < N; i++)
  {
    if (spamm_list_get(list, i) != i)
    {
      result = -1;
      printf("failed to sort, list[%u] = %i, should be %i...\n", i, spamm_list_get(list, i), i);
      break;
    }
  }

  return result;
}
