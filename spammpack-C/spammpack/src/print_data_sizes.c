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
  printf("&data at                     %p\n", &data);
  printf("&data.norm at                %p, offset = %lu\n", &data.norm, (unsigned long int) &data.norm - (unsigned long int) &data);
  printf("&data.block_dense at         %p, offset = %lu\n", &data.block_dense, (unsigned long int) &data.block_dense - (unsigned long int) &data);
  printf("&data.block_dense_dilated at %p, offset = %lu\n", &data.block_dense_dilated, (unsigned long int) &data.block_dense_dilated - (unsigned long int) &data);
}
