#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int i, j, k;
  struct spamm_t A, B, C;

  spamm_new(10, 10, 1, 1, 2, 2, 1e-10, &A);
  spamm_new(10, 10, 1, 1, 2, 2, 1e-10, &B);
  spamm_new(10, 10, 1, 1, 2, 2, 1e-10, &C);

  for (i = 0; i < 10; ++i) {
    for (j = 0; j < 10; ++j)
    {
      spamm_set(i, j, rand()/(double) RAND_MAX, &A);
      spamm_set(i, j, rand()/(double) RAND_MAX, &B);
    }
  }

  spamm_multiply(1.0, &A, &B, 1.0, &C);

  printf("A =\n");
  spamm_print_spamm(&A);
  printf("B =\n");
  spamm_print_spamm(&B);
  printf("C =\n");
  spamm_print_spamm(&C);

  return 0;
}
