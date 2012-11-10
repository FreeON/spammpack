#include <spamm.h>

#include <stdio.h>

int
main (int argc, char **argv)
{
  int result = 0;

  int dim;

  unsigned int i_test;

  unsigned int number_dimensions;
  unsigned int *i;
  unsigned int *N;

  unsigned int offset;

  for (number_dimensions = 1; number_dimensions <= 3; number_dimensions++)
  {
    i = calloc(number_dimensions, sizeof(unsigned int));
    N = calloc(number_dimensions, sizeof(unsigned int));

    for (i_test = 0; i_test < 100; i_test++)
    {
      for (dim = 0; dim < number_dimensions; dim++)
      {
        N[dim] = (unsigned int) (rand()/(double) RAND_MAX*100);
        i[dim] = (unsigned int) (rand()/(double) RAND_MAX*N[dim]);
      }

      switch (number_dimensions)
      {
        case 1:
          offset = i[0];
          if (spamm_index_row_major_3(number_dimensions, N, i) != offset)
          {
            for (dim = 0; dim < number_dimensions; dim++)
            {
              printf("N[%u] = %u, i[%u] = %u\n", dim, N[dim], dim, i[dim]);
            }
            SPAMM_FATAL("wrong index: offset = %u, found %u\n", offset,
                spamm_index_row_major_3(number_dimensions, N, i));
          }
          break;

        case 2:
          offset = i[1]+N[1]*i[0];
          if (spamm_index_row_major_3(number_dimensions, N, i) != offset)
          {
            for (dim = 0; dim < number_dimensions; dim++)
            {
              printf("N[%u] = %u, i[%u] = %u\n", dim, N[dim], dim, i[dim]);
            }
            SPAMM_FATAL("wrong index: offset = %u, found %u\n", offset,
                spamm_index_row_major_3(number_dimensions, N, i));
          }
          break;

        case 3:
          offset = i[2]+N[2]*(i[1]+N[1]*i[0]);
          if (spamm_index_row_major_3(number_dimensions, N, i) != offset)
          {
            for (dim = 0; dim < number_dimensions; dim++)
            {
              printf("N[%u] = %u, i[%u] = %u\n", dim, N[dim], dim, i[dim]);
            }
            SPAMM_FATAL("wrong index: offset = %u, found %u\n", offset,
                spamm_index_row_major_3(number_dimensions, N, i));
          }
          break;

        default:
          SPAMM_FATAL("FIXME\n");
          break;
      }
    }

    free(i);
    free(N);
  }

  return result;
}
