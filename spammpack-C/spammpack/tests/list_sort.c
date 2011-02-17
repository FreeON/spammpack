#include "spamm.h"
#include <math.h>
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
  float *node_norm;

  list = spamm_list_new(N);
  node_norm = calloc(sizeof(float), N);
  for (i = 0; i < N; i++)
  {
    spamm_list_set(list, i, rand());
    node_norm[i] = spamm_list_get(list, i);
  }

  spamm_list_sort(list, node_norm, spamm_list_compare_int);

  last_i = 0;
  last_value = spamm_list_get(list, 0);
  for (i = 0; i < N; i++)
  {
    value = spamm_list_get(list, i);
    if (fabs(value-node_norm[i]) > 1e-10)
    {
      printf("node norm was not sorted along with index\n");
      result = -1;
      break;
    }

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
  free(node_norm);

  return result;
}
