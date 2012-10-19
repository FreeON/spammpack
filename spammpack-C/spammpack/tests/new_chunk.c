#include <spamm.h>

int
main ()
{
  unsigned int number_dimensions = 2;
  unsigned int N_contiguous = 128;

  uint32_t *number_dimensions_pointer;
  uint32_t *N_contiguous_pointer;
  spamm_chunk_t *chunk;

  chunk = spamm_new_chunk(number_dimensions, N_contiguous);

  number_dimensions_pointer = spamm_chunk_get_number_dimensions(chunk);

  if (number_dimensions != *number_dimensions_pointer)
  {
    SPAMM_FATAL("number_dimensions mismatch\n");
  }

  N_contiguous_pointer = spamm_chunk_get_N_contiguous(chunk);

  if (N_contiguous != *N_contiguous_pointer)
  {
    SPAMM_FATAL("N_contiguous mismatch\n");
  }

  spamm_delete_chunk(&chunk);
}
