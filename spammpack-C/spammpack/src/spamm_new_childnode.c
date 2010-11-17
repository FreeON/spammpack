#include "spamm.h"
#include <math.h>

/** Allocate new child node.
 *
 * @param tier The tier this child node is on.
 * @param tree_depth The depth of the tree.
 *
 * @return A newly alloated child node.
 */
struct spamm_node_t *
spamm_new_childnode (const unsigned int tier,
    const unsigned int tree_depth,
    const unsigned int M_lower, const unsigned int M_upper,
    const unsigned int N_lower, const unsigned int N_upper,
    const unsigned int M_lower_kernel_tier, const unsigned int M_upper_kernel_tier,
    const unsigned int N_lower_kernel_tier, const unsigned int N_upper_kernel_tier,
    const unsigned int kernel_tier,
    floating_point_t *block_dense, floating_point_t *block_dense_dilated)
{
  unsigned int i, j;
  struct spamm_node_t *childnode = spamm_new_node();

  childnode->tier = tier;
  childnode->tree_depth = tree_depth;

  childnode->M_lower = M_lower;
  childnode->M_upper = M_upper;
  childnode->N_lower = N_lower;
  childnode->N_upper = N_upper;

  childnode->kernel_tier = kernel_tier;

  /* Check if we are at the kernel level. */
  if (childnode->tier == childnode->kernel_tier)
  {
    /* Store bounding box for kernel tier. */
    childnode->M_lower_kernel_tier = M_lower;
    childnode->M_upper_kernel_tier = M_upper;
    childnode->N_lower_kernel_tier = N_lower;
    childnode->N_upper_kernel_tier = N_upper;

    /* Reset newly allocated block to zero. */
    childnode->block_dense = (floating_point_t*) spamm_allocate(sizeof(floating_point_t)*SPAMM_N_KERNEL*SPAMM_N_KERNEL);
    for (i = 0; i < SPAMM_N_KERNEL; i++) {
      for (j = 0; j < SPAMM_N_KERNEL; j++)
      {
        childnode->block_dense[spamm_dense_index(i, j, SPAMM_N_KERNEL, SPAMM_N_KERNEL)] = 0.0;
      }
    }

    /* Allocate contiguous matrix block for dilated block. */
    childnode->block_dense_dilated = (floating_point_t*) spamm_allocate(sizeof(floating_point_t)*SPAMM_N_KERNEL*SPAMM_N_KERNEL*4);
    for (i = 0; i < SPAMM_N_KERNEL; i++) {
      for (j = 0; j < SPAMM_N_KERNEL; j++)
      {
        childnode->block_dense_dilated[spamm_dense_index(i, j, SPAMM_N_KERNEL, SPAMM_N_KERNEL)*4+0] = 0.0;
        childnode->block_dense_dilated[spamm_dense_index(i, j, SPAMM_N_KERNEL, SPAMM_N_KERNEL)*4+1] = 0.0;
        childnode->block_dense_dilated[spamm_dense_index(i, j, SPAMM_N_KERNEL, SPAMM_N_KERNEL)*4+2] = 0.0;
        childnode->block_dense_dilated[spamm_dense_index(i, j, SPAMM_N_KERNEL, SPAMM_N_KERNEL)*4+3] = 0.0;
      }
    }
  }

  else if (childnode->tier > childnode->kernel_tier)
  {
    childnode->M_lower_kernel_tier = M_lower_kernel_tier;
    childnode->M_upper_kernel_tier = M_upper_kernel_tier;
    childnode->N_lower_kernel_tier = N_lower_kernel_tier;
    childnode->N_upper_kernel_tier = N_upper_kernel_tier;

    if (childnode->tier < childnode->tree_depth)
    {
      /* We simply store the pointers to the dense blocks. */
      childnode->block_dense = block_dense;
      childnode->block_dense_dilated = block_dense_dilated;
    }

    else
    {
      /* Point into the contiguous matrix block. We use row-major order for
       * the blocks. */
      i = (M_lower-M_lower_kernel_tier)/SPAMM_N_BLOCK;
      j = (N_lower-N_lower_kernel_tier)/SPAMM_N_BLOCK;

      childnode->block_dense = block_dense+(i*SPAMM_N_KERNEL_BLOCK+j)*SPAMM_N_BLOCK*SPAMM_N_BLOCK;
      childnode->block_dense_dilated = block_dense_dilated+(i*SPAMM_N_KERNEL_BLOCK+j)*SPAMM_N_BLOCK*SPAMM_N_BLOCK*4;
    }
  }

  if (childnode->tier == childnode->kernel_tier)
  {
    /* Calculate linear Morton-ordered index of this node. We use the lower
     * corner of the index box for this calculation. The matrix indices
     * therefore jump by SPAMM_N_BLOCK.
     */
    spamm_print_node(childnode);
    childnode->index = spamm_linear_index_2D(childnode->M_lower, childnode->N_lower);
  }

  return childnode;
}
