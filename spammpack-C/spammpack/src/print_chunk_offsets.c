#include "spamm.h"

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

int
main (int argc, char **argv)
{
  unsigned int number_dimensions = 2;
  unsigned int N_contiguous = 128;
  unsigned int N_block = 4;

  unsigned int *N_lower;
  unsigned int *N_upper;

  spamm_chunk_t *chunk;

  unsigned int *N_lower_pointer;
  unsigned int *N_upper_pointer;
  float *A_pointer;
  float *A_dilated_pointer;
  float *norm_pointer;
  float *norm2_pointer;

  int dim;

  int c;
  const char *short_options = "hd:c:";
  const struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "dim", required_argument, NULL, 'd' },
    { "N_cont", required_argument, NULL, 'c' },
    { NULL, 0, NULL, 0 }
  };

  while (1)
  {
    c = getopt_long(argc, argv, short_options, long_options, NULL);

    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'h':
        printf("Usage:\n");
        printf("\n");
        printf("{ -h | --help }       This help\n");
        printf("{ -d | --dim } N      Set number of dimensions to N\n");
        printf("{ -c | --N_cont } N   Set N_contiguous to N\n");
        exit(0);
        break;

      case 'd':
        number_dimensions = strtol(optarg, NULL, 10);
        break;

      case 'c':
        N_contiguous = strtol(optarg, NULL, 10);
        break;

      default:
        printf("illegal option\n");
        exit(1);
        break;
    }
  }

  N_lower = calloc(number_dimensions, sizeof(unsigned int));
  N_upper = calloc(number_dimensions, sizeof(unsigned int));

  for (dim = 0; dim < number_dimensions; dim++)
  {
    N_lower[dim] = 0;
    N_upper[dim] = N_contiguous;
  }
  chunk = spamm_new_chunk(number_dimensions, N_block, N_lower, N_upper);

  free(N_lower);
  free(N_upper);

  printf("number_dimensions = %u\n", number_dimensions);
  printf("N_contiguous      = %u\n", N_contiguous);
  printf("sizeof(chunk)     = 0x%lx bytes\n",
      spamm_chunk_get_size(number_dimensions, N_block, N_lower, N_upper,
        &N_lower_pointer, &N_upper_pointer, &A_pointer, &A_dilated_pointer,
        &norm_pointer, &norm2_pointer));
  printf("\n");
  printf("&chunk->number_dimensions = 0x%lx\n", (void*) spamm_chunk_get_number_dimensions(chunk) - chunk);
  printf("&chunk->N_block           = 0x%lx\n", (void*) spamm_chunk_get_N_block(chunk) - chunk);
  printf("&chunk->N_lower           = 0x%lx\n", (void*) spamm_chunk_get_N_lower(chunk) - chunk);
  printf("&chunk->N_upper           = 0x%lx\n", (void*) spamm_chunk_get_N_upper(chunk) - chunk);
  printf("&chunk->A                 = 0x%lx\n", (void*) spamm_chunk_get_matrix(chunk) - chunk);
  printf("&chunk->A_dilated         = 0x%lx\n", (void*) spamm_chunk_get_matrix_dilated(chunk) - chunk);
}
