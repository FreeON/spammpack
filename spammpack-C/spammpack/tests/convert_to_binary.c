#include <stdio.h>
#include <stdlib.h>

  int
main (int argc, char **argv)
{
  double *A;
  double Aij;
  FILE *fd_in;
  FILE *fd_out;
  int N;
  int i, j;

  fd_in = fopen(argv[1], "r");
  N = 0;
  while (1)
  {
    if (fscanf(fd_in, "%i %i %le\n", &i, &j, &Aij) == EOF)
    {
      break;
    }
    if (i > N) N = i;
    if (j > N) N = j;
  }
  rewind(fd_in);

  printf("N = %i\n", N);
  A = calloc(N*N, sizeof(double));

  while (1)
  {
    if (fscanf(fd_in, "%i %i %le\n", &i, &j, &Aij) == EOF)
    {
      break;
    }
    A[(i-1)*N+(j-1)] = Aij;
  }
  fclose(fd_in);

  fd_out = fopen(argv[2], "wb");
  Aij = N;
  if (fwrite(&Aij, sizeof(double), 1, fd_out) < 1)
  {
    printf("error writing N to binary file\n");
    exit(1);
  }
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      if (fwrite(&A[i*N+j], sizeof(double), 1, fd_out) < 1)
      {
        printf("error writing matrix element to binary file, A(%i, %i)\n", i, j);
        exit(1);
      }
    }
  }
  fclose(fd_out);
}
