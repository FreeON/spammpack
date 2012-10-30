#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  const unsigned int N = 4;

  unsigned int dim;
  unsigned int number_dimensions = 3;
  unsigned int n, n_temp;
  unsigned int *i;

  i = calloc(number_dimensions, sizeof(unsigned int));

  for (n = 0; n < ipow(N, number_dimensions); n++)
  {
    n_temp = n;
    printf("i [");
    for (dim = 0; dim < number_dimensions; dim++)
    {
      i[dim] = n_temp%N;
      printf(" %u", i[dim]);
      if (dim+1 < number_dimensions)
      {
        printf(",");
      }
      n_temp /= N;
    }
    printf(" ] = %u\n", spamm_index_linear(number_dimensions, i));
  }

  return 0;
}
