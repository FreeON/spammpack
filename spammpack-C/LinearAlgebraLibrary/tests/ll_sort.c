#include <spamm.h>
#include <stdlib.h>

int
compare (const void *a_arg, const void *b_arg)
{
  int *a = (int*) a_arg;
  int *b = (int*) b_arg;

  if (*a < *b) { return -1; }
  else if (*a == *b) { return 0; }
  else { return +1; }
}

char *
data_to_string (const void *data)
{
  char *result;

  result = (char*) malloc(sizeof(char)*1000);
  snprintf(result, 1000, "%i", *((int*) data));
  return result;
}

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
    *data = 99-i;
    spamm_ll_append(data, list);
  }

  //spamm_ll_print(data_to_string, list);

  /* Sort. */
  spamm_ll_sort(compare, list);

  //spamm_ll_print(data_to_string, list);

  for (i = 0; i < 100; ++i)
  {
    data = spamm_ll_get(i, list);
    if (*data != i)
    {
      LOG("element %i should be %i but is %i\n", i, i, *data);
      return -1;
    }
  }

  /* Free memory. */
  spamm_ll_delete(free, &list);

  return 0;
}
