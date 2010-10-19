#include "spamm.h"
#include <stdio.h>
#include <stdlib.h>

/** \private Convert spamm_mm_t to string.
 *
 * @param data The data to convert.
 *
 * @return A string representing the data.
 */
char *
spamm_mm_print_data_to_string (const void *data)
{
  char *result = (char*) malloc(sizeof(char)*1000);

  snprintf(result, 1000, "%p", data);
  return result;
}

/** Print details of memory.
 *
 * @param memory The memory to print information about
 */
void
spamm_mm_print (const struct spamm_mm_t *memory)
{
  unsigned int i;
  struct spamm_ll_iterator_t *iterator;
  struct spamm_ll_node_t *node;
  struct spamm_mm_chunk_t *chunk;

  printf("details on memory: ");
  printf("%u chunks allocated\n", memory->chunks->number_elements);
  iterator = spamm_ll_iterator_new(memory->chunks);
  for (i = 0, node = spamm_ll_iterator_first(iterator); node != NULL; node = spamm_ll_iterator_next(iterator))
  {
    chunk = node->data;
    printf("chunk %u: chunksize = %u bytes, data from %p to %p, ", i++, chunk->chunksize, chunk->data, ((char*) chunk->data)+chunk->chunksize);
    printf("allocated = %u bytes\n", chunk->bytes_allocated);
    printf("allocated_start ");
    spamm_ll_print(NULL, chunk->allocated_start);
    printf("allocated_end   ");
    spamm_ll_print(NULL, chunk->allocated_end);
  }
}
