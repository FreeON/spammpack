#include "spamm.h"
#include <stdlib.h>

/** Swap spamm_multiply_stream_element_t.
 *
 * @param data1 The first element of type spamm_multiply_stream_element_t.
 * @param data2 The second element of type spamm_multiply_stream_element_t.
 */
void
spamm_swap_multiply_stream (void *data1, void *data2)
{
  struct spamm_multiply_stream_element_t *linear1 = data1;
  struct spamm_multiply_stream_element_t *linear2 = data2;

  if (linear1->M_A != linear2->M_A)
  {
    LOG2_FATAL("block dimensions do not match in M\n");
    exit(1);
  }

  if (linear1->N_A != linear2->N_A)
  {
    LOG2_FATAL("block dimensions do not match in N\n");
    exit(1);
  }

  if (linear1->M_B != linear2->M_B)
  {
    LOG2_FATAL("block dimensions do not match in M\n");
    exit(1);
  }

  if (linear1->N_B != linear2->N_B)
  {
    LOG2_FATAL("block dimensions do not match in N\n");
    exit(1);
  }

  if (linear1->M_C != linear2->M_C)
  {
    LOG2_FATAL("block dimensions do not match in M\n");
    exit(1);
  }

  if (linear1->N_C != linear2->N_C)
  {
    LOG2_FATAL("block dimensions do not match in N\n");
    exit(1);
  }

  spamm_swap_floating_point_t(&linear1->alpha, &linear2->alpha);
  spamm_swap_floating_point_t(&linear1->beta, &linear2->beta);

  spamm_swap_unsigned_int(&linear1->A_index, &linear2->A_index);
  spamm_swap_unsigned_int(&linear1->B_index, &linear2->B_index);
  spamm_swap_unsigned_int(&linear1->C_index, &linear2->C_index);

  spamm_swap_unsigned_int(&linear1->M_A, &linear2->M_A);
  spamm_swap_unsigned_int(&linear1->N_A, &linear2->N_A);
  spamm_swap_unsigned_int(&linear1->M_B, &linear2->M_B);
  spamm_swap_unsigned_int(&linear1->N_B, &linear2->N_B);
  spamm_swap_unsigned_int(&linear1->M_C, &linear2->M_C);
  spamm_swap_unsigned_int(&linear1->N_C, &linear2->N_C);

  spamm_swap_short_unsigned_int(&linear1->block_A_loaded_in_GPU, &linear2->block_A_loaded_in_GPU);
  spamm_swap_short_unsigned_int(&linear1->block_B_loaded_in_GPU, &linear2->block_B_loaded_in_GPU);
  spamm_swap_short_unsigned_int(&linear1->block_C_loaded_in_GPU, &linear2->block_C_loaded_in_GPU);

  spamm_swap_floating_point_t_pointer(&linear1->A_block_dense, &linear2->A_block_dense);
  spamm_swap_floating_point_t_pointer(&linear1->B_block_dense, &linear2->B_block_dense);
  spamm_swap_floating_point_t_pointer(&linear1->C_block_dense, &linear2->C_block_dense);

  spamm_swap_spamm_node_t_pointer(&linear1->C_node, &linear2->C_node);
}
