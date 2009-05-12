#include "lal.h"

int
lal_allocate_ (const int M, const int N, lal_matrix_t **A)
{
  return lal_allocate (M, N, A);
}
