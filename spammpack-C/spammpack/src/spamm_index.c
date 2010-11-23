#include "spamm.h"

unsigned int
spamm_dense_index (const unsigned int i, const unsigned int j)
{
  if (i >= SPAMM_N_KERNEL || i >= SPAMM_N_KERNEL)
  {
    fprintf(stderr, "i or j out of bounds\n");
    exit(1);
  }

  return i*SPAMM_N_KERNEL+j;
}
