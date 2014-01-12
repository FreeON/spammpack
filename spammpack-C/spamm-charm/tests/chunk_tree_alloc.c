#include "chunk_tree.h"

#include <stdio.h>

int
main ()
{
  const int N_chunk = 512;
  const int N_basic = 128;
  const int N = N_chunk;

  struct chunk_tree_t *chunk = chunk_tree_alloc(N_chunk, N_basic, N, 0, 0);

  printf("allocated %ld bytes\n", chunk_tree_sizeof(N_chunk, N_basic));

  chunk_tree_print(chunk, "chunk:\n");

  return 0;
}
