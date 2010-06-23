#include <spamm.h>
#include <stdlib.h>

int
main (int argc, char **argv)
{
  int i;
  int *data;
  struct spamm_ll_t *list;

  list = spamm_ll_new();

  for (i = 0; i < 100; ++i)
  {
    data = (int*) malloc(sizeof(int));
    *data = i;
    spamm_ll_append(data, list);
  }

  /* Delete the data. */
  while (list->number_elements > 0)
  {
    i = (int) (list->number_elements-1)*rand()/(double) RAND_MAX;
    data = spamm_ll_get(i, list);
    spamm_ll_delete_node(data, list);
  }
}
