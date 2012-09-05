#include "spamm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define N 100000
#define N_norm 100

int
main ()
{
  int result = 0;

  unsigned int index;
  float norm;

  unsigned int i, j, j_max;
  unsigned int last_i;
  unsigned int last_index;
  float last_norm;
  struct spamm_list_t *list;

  list = spamm_list_new(N);
  for (i = 0; i < N; i++)
  {
    index = rand();
    spamm_list_set(list, i, index, index);
  }

  /* Sort on index. */
  spamm_list_sort_index(list, spamm_list_compare_int);

  last_i = 0;
  last_index = spamm_list_get_index(list, 0);
  for (i = 0; i < N; i++)
  {
    index = spamm_list_get_index(list, i);
    norm = spamm_list_get_norm(list, i);
    if (fabs(index-norm) > 1e-10)
    {
      printf("node norm was not sorted along with index\n");
      result = -1;
      break;
    }

    if (index > last_index)
    {
      last_index = index;
      last_i = i;
    }

    else if (index < last_index)
    {
      result = -1;
      printf("failed to sort, list[%u] = %i, list[%u] = %i...\n", last_i, last_index, i, spamm_list_get_index(list, i));
      break;
    }
  }

  /* Sort on norm. */
  for (j = 0; j < N; j += N_norm)
  {
    j_max = (j+N_norm <= N ? j+N_norm : N);
    spamm_list_sort_norm(list, j, j_max);

    last_i = j;
    last_norm = spamm_list_get_norm(list, j);
    for (i = j; i < j_max; i++)
    {
      index = spamm_list_get_index(list, i);
      norm = spamm_list_get_norm(list, i);
      if (fabs(index-norm) > 1e-10)
      {
        printf("node norm was not sorted along with index\n");
        result = -1;
        break;
      }

      if (norm < last_norm)
      {
        last_norm = norm;
        last_i = i;
      }

      else if (norm > last_norm)
      {
        result = -1;
        printf("[%u, %u[: failed to sort, list[%u] = %f, list[%u] = %f...\n",
            j, j_max, last_i, last_norm, i, spamm_list_get_norm(list, i));
        break;
      }
    }
  }

  spamm_list_delete(&list);

  return result;
}
