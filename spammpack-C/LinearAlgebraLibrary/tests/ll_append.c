#include <spamm.h>
#include <stdlib.h>

int
main ()
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

  if (list->number_elements != 100)
  {
    LOG("wrong number of elements, should be 100, is reported to be %u\n", list->number_elements);
    return -1;
  }

  for (i = 0; i < 100; ++i)
  {
    data = spamm_ll_get(i, list);
    if (i != *data)
    {
      LOG("wrong element %i != %i\n", i, *data);
      return -1;
    }
  }

  /* Free memory for list. */
  spamm_ll_delete(list);

  return 0;
}
