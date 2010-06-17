#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int N = 10;
  int chunksize = 12;

  int i, j;
  int **i_pointer;
  struct spamm_mm_t *memory;
  struct spamm_ll_node_t *node;
  struct spamm_ll_iterator_t *iterator;
  struct spamm_mm_chunk_t *chunk;

  /* Allocate new memory. */
  memory = spamm_mm_new(chunksize);

  /* Store. */
  i_pointer = (int**) malloc(sizeof(int*)*N);
  for (i = 0; i < N; ++i)
  {
    i_pointer[i] = spamm_mm_allocate(sizeof(int), memory);
    *(i_pointer[i]) = i;
  }

  /* Retrieve. */
  iterator = spamm_ll_iterator_new(memory->chunks);
  for (i = 0, node = spamm_ll_iterator_first(iterator); node != NULL; node = spamm_ll_iterator_next(iterator))
  {
    chunk = node->data;
    for (j = 0; j < chunk->allocated_start->number_elements; ++j)
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
