#include "spamm.h"
#include <stdio.h>

int
main ()
{
  struct spamm_data_t data;
  struct spamm_multiply_stream_t stream_element;

  printf("sizeof(struct spamm_multiply_stream_t) = %lu\n", sizeof(struct spamm_multiply_stream_t));
  printf("&stream_element at   %p\n", &stream_element);
  printf("&stream_element.A at %p, offset = %lu\n", &stream_element.A, (unsigned long int) &stream_element.A - (unsigned long int) &stream_element);
  printf("&stream_element.B at %p, offset = %lu\n", &stream_element.B, (unsigned long int) &stream_element.B - (unsigned long int) &stream_element);
  printf("&stream_element.C at %p, offset = %lu\n", &stream_element.C, (unsigned long int) &stream_element.C - (unsigned long int) &stream_element);

  printf("sizeof(struct spamm_data_t) = %lu\n", sizeof(struct spamm_data_t));
  printf("&data at                       %p\n", &data);
  printf("&data.norm at                  %p, offset_norm = %lu\n", &data.norm, (unsigned long int) &data.norm - (unsigned long int) &data);
#ifdef SPAMM_KERNEL_IMPLEMENTATION_Z_CURVE_SSE4_1
  printf("&data.norm_upper at            %p, offset_norm_upper = %lu\n", &data.norm_upper, (unsigned long int) &data.norm_upper - (unsigned long int) &data);
  printf("&data.norm_upper_transpose at  %p, offset_norm_upper_transpose = %lu\n", &data.norm_upper_transpose, (unsigned long int) &data.norm_upper_transpose - (unsigned long int) &data);
#endif
  printf("&data.block_dense at           %p, offset_block_dense = %lu\n", &data.block_dense, (unsigned long int) &data.block_dense - (unsigned long int) &data);
  printf("&data.block_dense_dilated at   %p, offset_block_dense_dilated = %lu\n", &data.block_dense_dilated, (unsigned long int) &data.block_dense_dilated - (unsigned long int) &data);
#ifdef SPAMM_USE_TRANSPOSE
  printf("&data.block_dense_transpose at %p, offset_block_dense_transpose = %lu\n", &data.block_dense_transpose, (unsigned long int) &data.block_dense_transpose - (unsigned long int) &data);
#endif
}
