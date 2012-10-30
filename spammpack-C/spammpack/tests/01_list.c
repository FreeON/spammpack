#include "spamm_list.h"
#include <stdio.h>
#include <stdlib.h>

#define N 1000

int
main ()
{
  unsigned int i;
  struct spamm_list_t *list;

  list = spamm_list_new(N);

  if (spamm_list_length(list) != N)
  {
    printf("list length incorrect\n");
    exit(1);
  }

  for (i = 0; i < N; i++)
  {
    spamm_list_set(list, i, i, i);
  }

  for (i = 0; i < N; i++)
  {
    if (spamm_list_get_index(list, i) != i)
    {
      printf("incorrect list entry for i = %u, found %u\n", i, spamm_list_get_index(list, i));
      exit(1);
    }
  }

  spamm_list_delete(&list);

  return 0;
}
