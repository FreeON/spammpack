#include "spamm.h"

unsigned int
spamm_index_row_major (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N)
{
  return i*N+j;
}

unsigned int
spamm_index_kernel_block (const unsigned int i, const unsigned int j)
{
  unsigned int offset;

  offset = SPAMM_N_BLOCK*SPAMM_N_BLOCK*spamm_index_row_major(i/SPAMM_N_BLOCK,
      j/SPAMM_N_BLOCK, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)
    + spamm_index_row_major(i%SPAMM_N_BLOCK, j%SPAMM_N_BLOCK, SPAMM_N_BLOCK, SPAMM_N_BLOCK);

  return offset;
}
