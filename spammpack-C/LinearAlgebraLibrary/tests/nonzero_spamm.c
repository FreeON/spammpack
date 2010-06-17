#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int max_size = 5;
  int M[5] = { 2, 10, 15, 200, 300 };
  int N[5] = { 2, 10,  8, 200, 100 };

  int max_block = 11;
  int M_block[11] = { 1, 1, 1, 2, 2, 2, 5, 8, 8, 10, 100 };
  int N_block[11] = { 1, 2, 3, 1, 2, 5, 2, 8, 2, 10, 100 };

  int max_child = 5;
  int M_child[5] = { 2, 2, 2, 3, 4 };
  int N_child[5] = { 2, 3, 4, 3, 4 };

  int max_fill = 4;
  double fill[4] = { 0.01, 0.2, 0.5, 1.0 };

  int i, j, k;
  int i_size;
  int i_block;
  int i_child;
  int i_fill;

  int number_nonzero;
  struct spamm_t A;

  int result = 0;

  for (i_fill = 0; i_fill < max_fill; ++i_fill) {
    for (i_size = 0; i_size < max_size; ++i_size) {
      for (i_block = 0; i_block < max_block; ++i_block) {
        for (i_child = 0; i_child < max_child; ++i_child)
        {
          spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &A);
          number_nonzero = (int) (M[i_size]*N[i_size]*fill[i_fill]);

          //printf("[nonzero_spamm] %ix%i matrix, %.1f%% full, padded %ix%i, i_block dimensions %ix%i, i_child dimensions %ix%i, %i nonzeros\n",
          //    M[i_size], N[i_size], fill[i_fill]*100, A.M_padded, A.N_padded, M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child],
          //    number_nonzero);

          for (k = 0; k < number_nonzero; ++k)
          {
            while (1)
            {
              i = (int) (rand()/((double) RAND_MAX+1)*M[i_size]);
              j = (int) (rand()/((double) RAND_MAX+1)*N[i_size]);
              if (spamm_get(i, j, &A) == 0.0)
              {
                spamm_set(i, j, 1.0, &A);
                break;
              }
            }
          }

          if (spamm_number_nonzero(&A) != number_nonzero)
          {
            spamm_log("set %i nonzeros, but found %i\n", __FILE__, __LINE__, number_nonzero, spamm_number_nonzero(&A));
            result = 1;
          }
        }
      }
    }
  }

  spamm_delete(&A);

  return result;
}
