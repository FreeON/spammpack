#include <spamm.h>
#include <spamm_ll.h>
#include <spamm_mm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int N = 10;
  int chunksize = 12;

  int i, j;
  int **i_pointer;
  struct spamm_ll_t *memory;
  struct spamm_ll_node_t *node;
  struct spamm_mm_t *chunk;

  memory = spamm_mm_initialize(chunksize);

  /* Store. */
  i_pointer = (int**) malloc(sizeof(int*)*N);
  for (i = 0; i < 10; ++i)
  {
    i_pointer[i] = spamm_mm_allocate(sizeof(int), memory);
    *(i_pointer[i]) = i;
  }

  /* Retrieve. */
  for (i = 0, node = memory->first; node != NULL; node = node->next)
  {
    chunk = node->data;
    for (j = 0; j < chunk->allocated_start.number_elements; ++j)
    {
      if (i != *(i_pointer[i]))
      {
        LOG("value mismatch: found %i but should be %i\n", *(i_pointer[i]), i);
        return -1;
      }
      i++;
    }
  }

  return 0;
}
