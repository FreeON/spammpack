#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** \private Pack a subtree.
 *
 * Packing means that the lowest tiers of the matrix tree are converted into a
 * linear quadtree and packed into chunks of contiguous memory. The matrix
 * data in the chunks is sorted according to the data's linear index.
 *
 * @param linear_tier The tier of linear quadtree storage.
 * @param chunksize The size in bytes of the data chunks.
 * @param mask The mask to apply to the linear index before sorting.
 * @param linear_tree The spamm_mm_t memory that holds the linear tree.
 * @param node The node to pack.
 */
void
spamm_tree_pack_subtree (const unsigned int linear_tier, const unsigned int chunksize,
    const enum spamm_linear_mask_t mask, struct spamm_ll_t *linear_tree,
    struct spamm_node_t *node)
{
  int i, j;
  struct spamm_linear_quadtree_t *linear_block;

  node->linear_tier = linear_tier;

  if (node->tier == linear_tier)
  {
    /* Pack. */
    LOG("reached tier %u, packing\n", node->tier);
    linear_tree = spamm_mm_initialize(chunksize);
    node->linear_quadtree = linear_tree;
    LOG("linear_tree at %p\n", linear_tree);
  }

  /* Recurse more. */
  if (node->child != NULL)
  {
    LOG("linear_tree at %p\n", linear_tree);
    for (i = 0; i < node->M_child; ++i) {
      for (j = 0; j < node->N_child; ++j)
      {
        spamm_tree_pack_subtree(linear_tier, chunksize, mask, linear_tree, node->child[spamm_dense_index(i, j, node->M_child, node->N_child)]);
      }
    }
  }

  else
  {
    LOG("linear_tree at %p\n", linear_tree);
    if (node->block_dense != NULL)
    {
      /* Copy data into linear tree. */
      LOG("tier %u, storing datablock with index %u\n", node->tier, node->index);
      linear_block = (struct spamm_linear_quadtree_t*) spamm_mm_allocate(sizeof(struct spamm_linear_quadtree_t)+node->M_block*node->N_block*sizeof(float_t)+8, linear_tree);

      LOG("linear_block at %p\n", linear_block);

      /* Move linear_block->block_dense pointer to just after the pointer
       * itself. Since we allocated enough space in spamm_mm_allocate() to fit
       * the pointer and the dense matrix block, we can do that without having
       * to allocate the dense data block separately.
       */
      linear_block->block_dense = (float_t*) (&(linear_block->block_dense)+1);
      linear_block->index = node->index;

      /* Copy data. */
      for (i = 0; i < node->M_block; ++i) {
        for (j = 0; j < node->N_block; ++j)
        {
          linear_block->block_dense[spamm_dense_index(i, j, node->M_block, node->N_block)] = node->block_dense[spamm_dense_index(i, j, node->M_block, node->N_block)];
        }
      }
    }
  }

  if (node->tier == linear_tier)
  {
    /* Delete tree, which is replaced now by linear quadtree. */
    spamm_log("done packing, deleting original tree\n", __FILE__, __LINE__);
    spamm_delete_node(node);
  }
}

/** Pack a matrix tree.
 *
 * Packing means that the lowest tiers of the matrix tree are converted into a
 * linear quadtree and packed into chunks of contiguous memory. The matrix
 * data in the chunks is sorted according to the data's linear index.
 *
 * @param linear_tier The tier of linear quadtree storage.
 * @param chunksize The size in bytes of the data chunks.
 * @param mask The mask to apply to the linear index before sorting.
 * @param A The matrix to pack.
 *
 * \bug Support for the linear quadtree is not implemented everywhere yet.
 *      Trees that are packed can not not be used in other functions.
 */
void
spamm_tree_pack (const unsigned int linear_tier, const unsigned int chunksize,
    const enum spamm_linear_mask_t mask, struct spamm_t *A)
{
  assert(A != NULL);

  if (linear_tier > A->tree_depth)
  {
    LOG("tree depth is only %u, while linear_tier is %u\n", A->tree_depth, linear_tier);
    exit(1);
  }
  A->linear_tier = linear_tier;

  spamm_log("packing tree\n", __FILE__, __LINE__);

  /* Recurse tree and convert to linear trees any subtrees below linear_tier.
   */
  if (A->root != NULL)
  {
    spamm_tree_pack_subtree(linear_tier, chunksize, mask, NULL, A->root);
  }

  else
  {
    spamm_log("A is empty, nothing to pack\n", __FILE__, __LINE__);
  }
}
