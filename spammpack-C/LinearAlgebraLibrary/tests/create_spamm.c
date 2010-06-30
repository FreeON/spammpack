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

  int i, j;
  int i_size;
  int i_block;
  int i_child;
  int i_fill;

  floating_point_t *A_dense;
  struct spamm_t A;
  struct spamm_tree_stats_t stats;

  int result = 0;

  for (i_fill = 0; i_fill < max_fill; ++i_fill)
  {
    //i_fill = 3;
    for (i_size = 0; i_size < max_size; ++i_size)
    {
      //i_size = 1;
      A_dense = (floating_point_t*) malloc(sizeof(floating_point_t)*M[i_size]*N[i_size]);
      for (i = 0; i < M[i_size]; ++i) {
        for (j = 0; j < N[i_size]; ++j)
        {
          A_dense[spamm_dense_index(i, j, M[i_size], N[i_size])] = (rand()/(floating_point_t) RAND_MAX > (1-fill[i_fill]) ? rand()/(floating_point_t) RAND_MAX : 0.0);
        }
      }

      //printf("A_dense =\n");
      //spamm_print_dense(M[i_size], N[i_size], A_dense);

      for (i_block = 0; i_block < max_block; ++i_block)
      {
        //i_block = 0;
        for (i_child = 0; i_child < max_child; ++i_child)
        {
          //i_child = 1;
          spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &A);
          //printf("[create_spamm] %ix%i matrix, %.1f%% full, padded %ix%i, i_block dimensions %ix%i, i_child dimensions %ix%i, ",
          //    M[i_size], N[i_size], fill[i_fill]*100, A.M_padded, A.N_padded, M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child]);
          spamm_dense_to_spamm(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, A_dense, &A);
          spamm_tree_stats(&stats, &A);
          //printf("nodes = %i, blocks = %i, tree = %i bytes, blocks = %i bytes (%1.1f%%)\n", stats.number_nodes, stats.number_dense_blocks,
          //    stats.memory_tree, stats.memory_dense_blocks, (stats.memory_tree+stats.memory_dense_blocks)/(double) (M[i_size]*N[i_size]*sizeof(double))*100);

          //printf("A tree =\n");
          //spamm_print_tree(&A);
          //printf("A =\n");
          //spamm_print_spamm(&A);

          //printf("i_fill = %i, i_size = %i, i_block = %i, i_child = %i\n", i_fill, i_size, i_block, i_child);

          for (i = 0; i < M[i_size]; ++i) {
            for (j = 0; j < N[i_size]; ++j)
            {
              if (A_dense[spamm_dense_index(i, j, M[i_size], N[i_size])] != spamm_get(i, j, &A))
              {
                printf("[create_spamm] mismatch: (A_dense[%i][%i] = %e) != (A[%i][%i] = %e)\n", i, j, A_dense[spamm_dense_index(i, j, M[i_size], N[i_size])], i, j, spamm_get(i, j, &A));
                result = 1;
              }
            }
          }

          spamm_delete(&A);
        }
      }

      free(A_dense);
    }
  }

  return result;
}
