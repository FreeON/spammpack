#include "spamm.h"

int
spamm_dense_index (const int i, const int j, const int stride)
{
  return i*stride+j;
}
