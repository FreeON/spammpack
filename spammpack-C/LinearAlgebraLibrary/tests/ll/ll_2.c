#include <spamm.h>
#include <math.h>
#include <stdlib.h>

int
main (int argc, char **argv)
{
  int N = 10;

  int i, j, temp;
  int *data;
  int *removal;
  struct spamm_ll_iterator_t *iterator;
  struct spamm_ll_node_t *node;
  struct spamm_ll_t *list;

  list = spamm_ll_new();

  //spamm_set_loglevel(debug);

  for (i = 0; i < N; ++i)
  {
    data = (int*) malloc(sizeof(int));
    *data = i;
    spamm_ll_append(data, list);
  }

  /* Create removal list. */
  removal = (int*) malloc(sizeof(int)*N);
  for (i = 0; i < N; ++i)
  {
    removal[i] = i;
  }

  /* Shuffle removal list. */
  for (i = 0; i < N-1; ++i) {
    for (j = i+1; j < N; ++j)
    {
      if (rand()/(double) RAND_MAX > 0.5)
      {
        temp = removal[i];
        removal[i] = removal[j];
        removal[j] = temp;
      }
    }
  }

  /* Delete the data. */
  for (i = 0; i < N; ++i)
  {
    LOG_DEBUG("deleting i = %i\n", removal[i]);

    iterator = spamm_ll_iterator_new(list);
    for (node = spamm_ll_iterator_first(iterator); node != NULL; node = spamm_ll_iterator_next(iterator))
    {
      data = node->data;
      if (*data == removal[i])
      {
        break;
      }
    }
    spamm_ll_iterator_delete(&iterator);

    if (node != NULL)
    {
      spamm_ll_delete_node(free, node, list);
      LOG_DEBUG("list has %u elements\n", list->number_elements);
    }

    else
    {
      LOG2_FATAL("node == NULL?\n");
      return -1;
    }
  }

  if (list->number_elements != 0)
  {
    LOG2_FATAL("list still has elements\n");
    return -1;
  }

  free(removal);
  spamm_ll_delete(free, &list);

  return 0;
}
