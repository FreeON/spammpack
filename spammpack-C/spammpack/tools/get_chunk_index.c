#include "spamm.h"

#include <stdio.h>

void
main ()
{
  int i[2];
  const unsigned int N_lower[2] = { 0, 0 };
  const unsigned int N_upper[2] = { 16, 16 };
  unsigned int offset;

  SPAMM_INFO("N_lower = { %u, %u }\n", N_lower[0], N_lower[1]);
  SPAMM_INFO("N_upper = { %u, %u }\n", N_upper[0], N_upper[1]);

  while(1)
  {
    printf("Enter i j: "); fflush(stdout);

    if(fscanf(stdin, "%i %i", &i[0], &i[1]) != 2)
    {
      break;
    }

    if(i[0] < 0 || i[1] < 0)
    {
      break;
    }

    offset = spamm_chunk_matrix_index(2, 1, N_lower, N_upper, i);
    SPAMM_INFO("i = %i, j = %i, chunk_offset = %u, # kernel blocks = %u\n", i[0],
        i[1], offset, offset/(spamm_chunk_get_kernel_size()*spamm_chunk_get_kernel_size()));
  }
}
