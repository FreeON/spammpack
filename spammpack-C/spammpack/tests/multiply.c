#include "spamm.h"

#include <stdio.h>
#include <stdlib.h>

int
main (int argc, char *argv[])
{
  int result = 0;
  const int M = 103;
  const int N = 129;

  unsigned int i, j;
  float Aij;

  FILE *fd;

  const char *filename_A = "testmatrix_1.coor";

  struct spamm_matrix_t *A;

  /* Allocate matrices. */
  A = spamm_new(M, N);

  /* Load matrix from file. */
  if ((fd = fopen(filename_A, "r")) == NULL)
  {
    printf("error opening matrix file \"%s\"\n", filename_A);
    exit(1);
  }

  while (fscanf(fd, "%u %u %e\n", &i, &j, &Aij) == 3)
  {
    spamm_set(i, j, Aij, A);
  }

  if (fclose(fd) != 0)
  {
    printf("error closing matrix file\n");
    exit(1);
  }

  return result;
}
