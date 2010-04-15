#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int max_size = 4;
  int M[4] = { 2, 10, 1000, 5000 };
  int N[4] = { 2, 10, 1000, 5000};

  int max_block = 11;
  int M_block[11] = { 1, 1, 1, 2, 2, 2, 5, 8, 8, 10, 100 };
  int N_block[11] = { 1, 2, 3, 1, 2, 5, 2, 8, 2, 10, 100 };

  int max_child = 3;
  int M_child[3] = { 2, 3, 4 };
  int N_child[3] = { 2, 3, 4 };

  int i, j, k, l, m;
  double *A_dense;
  struct spamm_t A;

  for (k = 0; k < max_size; ++k)
  {
    A_dense = (double*) malloc(sizeof(double)*M[k]*N[k]);
    for (i = 0; i < M[k]; ++i) {
      for (j = 0; j < N[k]; ++j)
      {
        A_dense[spamm_dense_index(i, j, N[k])] = rand()/(double) RAND_MAX;
      }
    }

    for (l = 0; l < max_block; ++l)
    {
      for (m = 0; m < max_child; ++m)
      {
        spamm_new(M[k], N[k], M_block[l], N_block[l], M_child[m], N_child[m], &A);
        printf("%ix%i matrix, padded %ix%i, block dimensions %ix%i, child dimensions %ix%i\n", M[k], N[k], A.M_padded, A.N_padded, M_block[l], N_block[l], M_child[m], N_child[m]);
        spamm_dense_to_spamm(M[k], N[k], M_block[l], N_block[l], M_child[m], N_child[m], A_dense, &A);

        for (i = 0; i < M[k]; ++i) {
          for (j = 0; j < N[k]; ++j)
          {
            if (A_dense[spamm_dense_index(i, j, N[k])] != spamm_get(i, j, &A))
            {
              printf("mismatch: (A_dense[%i][%i] = %e) != (A[%i][%i] = %e)\n", i, j, A_dense[spamm_dense_index(i, j, N[k])], i, j, spamm_get(i, j, &A));
              exit(1);
            }
          }
        }

        spamm_delete(&A);
      }
    }

    free(A_dense);
  }

  fprintf(stderr, "test ok\n");
  return 0;
}
