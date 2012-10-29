#include <spamm.h>

int
main ()
{
  unsigned int number_dimensions = 2;
  unsigned int N_contiguous = 128;
  unsigned int N_block = 4;

  unsigned int *N;
  unsigned int *N_lower;
  unsigned int *N_upper;

  unsigned int *number_dimensions_pointer;

  spamm_chunk_t *chunk;

  int dim;

  N = calloc(number_dimensions, sizeof(unsigned int));
  N_lower = calloc(number_dimensions, sizeof(unsigned int));
  N_upper = calloc(number_dimensions, sizeof(unsigned int));

  for (dim = 0; dim < number_dimensions; dim++)
  {
    N[dim] = N_contiguous;
    N_lower[dim] = 0;
    N_upper[dim] = N_contiguous;
  }
  chunk = spamm_new_chunk(number_dimensions, N_block, N, N_lower, N_upper);

  free(N_lower);
  free(N_upper);

  number_dimensions_pointer = spamm_chunk_get_number_dimensions(chunk);

  if (number_dimensions != *number_dimensions_pointer)
  {
    SPAMM_FATAL("number_dimensions mismatch\n");
  }

  if (N_contiguous != spamm_chunk_get_N_contiguous(chunk))
  {
    SPAMM_FATAL("N_contiguous mismatch\n");
  }

  spamm_delete_chunk(&chunk);
}
