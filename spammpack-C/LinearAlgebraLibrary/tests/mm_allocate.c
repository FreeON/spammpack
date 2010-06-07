#include <spamm_ll.h>
#include <spamm_mm.h>
#include <stdio.h>

int
main ()
{
  int i, j;
  int *i_pointer;
  struct spamm_ll_t *memory;
  struct spamm_ll_node_t *node;
  struct spamm_mm_t *chunk;

  memory = spamm_mm_initialize(12);

  /* Store. */
  for (i = 0; i < 10; ++i)
  {
    i_pointer = spamm_mm_allocate(sizeof(int), memory);
    *i_pointer = i;
  }

  /* Retrieve. */
  for (i = 0, node = memory->first; node != NULL; node = node->next)
  {
    chunk = node->data;
    i_pointer = chunk->data;
    for (j = 0; j < chunk->allocated_start.number_elements; ++j)
    {
      if (i != i_pointer[j])
      {
        spamm_log("value mismatch: found %i but should be %i\n", i_pointer[j], i);
        return -1;
      }
      //printf("%i = %i (%p)\n", i, i_pointer[j], &(i_pointer[j]));
      i++;
    }
  }

  return 0;
}
