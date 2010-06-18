#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int max_size = 5;
  int M[] = { 2, 10, 15, 200, 300 };
  int N[] = { 2, 10,  8, 200, 100 };

  int max_block = 10;
  int M_block[] = { 1, 1, 1, 2, 2, 2, 5, 8, 8, 10 };
  int N_block[] = { 1, 2, 3, 1, 2, 5, 2, 8, 2, 10 };

  int max_child = 1;
  int M_child[] = { 2 };
  int N_child[] = { 2 };

  int max_fill_A = 4;
  double fill_A[] = { 0.01, 0.2, 0.5, 1.0 };

  int max_fill_B = 4;
  double fill_B[] = { 0.01, 0.2, 0.5, 1.0 };

  int number_linear_tiers = 3;
  int linear_tier[] = { 0, 1, 3 };

  int number_chunks = 3;
  int chunksize[] = { 500, 1000, 10000 };

  int i, j;
  int i_size;
  int i_block;
  int i_child;
  int i_fill_A;
  int i_fill_B;
  int i_linear;
  int i_chunk;

  struct spamm_t A, B, B_test;
  struct spamm_tree_stats_t stats;

  int result = 0;

  //spamm_set_loglevel(debug);

  for (i_fill_A = 0; i_fill_A < max_fill_A; ++i_fill_A) {
    for (i_fill_B = 0; i_fill_B < max_fill_B; ++i_fill_B) {
      for (i_size = 0; i_size < max_size; ++i_size) {
        for (i_block = 0; i_block < max_block; ++i_block) {
          for (i_child = 0; i_child < max_child; ++i_child) {
            for (i_linear = 0; i_linear < number_linear_tiers; ++i_linear) {
              for (i_chunk = 0; i_chunk < number_chunks; ++i_chunk)
              {
                spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &A);
                spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &B);
                spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &B_test);

                //LOG2_INFO("starting test...\n");
                printf("%ix%i matrix, A = %.1f%% full, B = %.1f%% full, padded %ix%i, i_block dimensions %ix%i, i_child dimensions %ix%i\n",
                    M[i_size], N[i_size], fill_A[i_fill_A]*100,
                    fill_B[i_fill_B]*100, A.M_padded, A.N_padded,
                    M_block[i_block], N_block[i_block], M_child[i_child],
                    N_child[i_child]);

                for (i = 0; i < M[i_size]; ++i) {
                  for (j = 0; j < N[i_size]; ++j)
                  {
                    spamm_set(i, j, (rand()/(double) RAND_MAX > (1-fill_A[i_fill_A]) ? rand()/(double) RAND_MAX : 0.0), &A);
                    spamm_set(i, j, (rand()/(double) RAND_MAX > (1-fill_B[i_fill_B]) ? rand()/(double) RAND_MAX : 0.0), &B);
                  }
                }

                //printf("A =\n");
                //spamm_print_tree(&A);
                //spamm_tree_stats(&stats, &A);
                //LOG_INFO("nodes = %i, blocks = %i, tree = %i bytes, blocks = %i bytes (%1.1f%%)\n", stats.number_nodes, stats.number_dense_blocks,
                //    stats.memory_tree, stats.memory_dense_blocks,
                //    (stats.memory_tree+stats.memory_dense_blocks)/(double) (M[i_size]*N[i_size]*sizeof(double))*100);
                //fflush(stdout);

                //printf("B =\n");
                //spamm_print_tree(&B);
                //spamm_tree_stats(&stats, &B);
                //LOG_INFO("nodes = %i, blocks = %i, tree = %i bytes, blocks = %i bytes (%1.1f%%)\n", stats.number_nodes, stats.number_dense_blocks,
                //    stats.memory_tree, stats.memory_dense_blocks,
                //    (stats.memory_tree+stats.memory_dense_blocks)/(double) (M[i_size]*N[i_size]*sizeof(double))*100);
                //fflush(stdout);

                /* Pack linear quadtrees. */
                spamm_tree_pack(linear_tier[i_linear], chunksize[i_chunk], i_mask, &A);
                spamm_tree_pack(linear_tier[i_linear], chunksize[i_chunk], i_mask, &B);

                //printf("A =\n");
                //spamm_print_tree(&A);
                //spamm_tree_stats(&stats, &A);
                //LOG_DEBUG("nodes = %i, blocks = %i, tree = %i bytes, blocks = %i bytes (%1.1f%%)\n", stats.number_nodes, stats.number_dense_blocks,
                //    stats.memory_tree, stats.memory_dense_blocks, (stats.memory_tree+stats.memory_dense_blocks)/(double) (M[i_size]*N[i_size]*sizeof(double))*100);
                //fflush(stdout);

                //printf("B =\n");
                //spamm_print_tree(&B);
                //spamm_tree_stats(&stats, &B);
                //LOG_DEBUG("nodes = %i, blocks = %i, tree = %i bytes, blocks = %i bytes (%1.1f%%)\n", stats.number_nodes, stats.number_dense_blocks,
                //    stats.memory_tree, stats.memory_dense_blocks, (stats.memory_tree+stats.memory_dense_blocks)/(double) (M[i_size]*N[i_size]*sizeof(double))*100);
                //fflush(stdout);

                //printf("A =\n");
                //spamm_print_spamm(&A);
                //printf("B =\n");
                //spamm_print_spamm(&B);

                for (i = 0; i < M[i_size]; ++i) {
                  for (j = 0; j < N[i_size]; ++j)
                  {
                    spamm_set(i, j, spamm_get(i, j, &A)+spamm_get(i, j, &B), &B_test);
                  }
                }

                //printf("B_test =\n");
                //spamm_print_spamm(&B_test);

                spamm_add(1.0, &A, 1.0, &B);

                //printf("B =\n");
                //spamm_print_spamm(&B);

                for (i = 0; i < M[i_size]; ++i) {
                  for (j = 0; j < N[i_size]; ++j)
                  {
                    if (spamm_get(i, j, &B) != spamm_get(i, j, &B_test))
                    {
                      LOG_FATAL("mismatch: (B_test[%i][%i] = %e) != (B[%i][%i] = %e)\n", i, j, spamm_get(i, j, &B_test), i, j, spamm_get(i, j, &B));
                      return -1;
                    }
                  }
                }
                LOG2_DEBUG("comparison ok\n");

                spamm_delete(&A);
                spamm_delete(&B);
                spamm_delete(&B_test);
              }
            }
          }
        }
      }
    }
  }

  return result;
}
