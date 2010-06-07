#include "spamm.h"
#include "spamm_mm.h"
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
    const enum spamm_linear_mask_t mask, struct spamm_mm_t *linear_tree,
    struct spamm_node_t *node)
{
  int i, j;
  struct spamm_linear_quadtree_t *linear_block;

  node->linear_tier = linear_tier;

  if (node->tier == linear_tier)
  {
    /* Pack. */
    LOG("reached tier %u, packing\n", node->tier);
    //linear_tree = spamm_mm_initialize(chunksize);
  }

  /* Recurse more. */
  if (node->child != NULL)
  {
    for (i = 0; i < node->M_child; ++i) {
      for (j = 0; j < node->N_child; ++j)
      {
        spamm_tree_pack_subtree(linear_tier, chunksize, mask, linear_tree, node->child[spamm_dense_index(i, j, node->M_child, node->N_child)]);
      }
    }
  }

  else
  {
    if (node->block_dense != NULL)
    {
      /* Copy data into linear tree. */
      //linear_block = (struct spamm_linear_quadtree_t*) spamm_mm_allocate(sizeof(struct spamm_linear_quadtree_t)+node->M_block*node->N_block*sizeof(float_t), linear_tree);
      LOG("storing datablock with index %u\n", node->index);
      linear_block = (struct spamm_linear_quadtree_t*) malloc(sizeof(struct spamm_linear_quadtree_t));
      linear_block->block_dense = (float_t*) malloc(sizeof(float_t)*node->M_block*node->N_block);
      linear_block->index = node->index;
      for (i = 0; i < node->M_block; ++i) {
        for (j = 0; j < node->N_block; ++j)
        {
          linear_block->block_dense[spamm_dense_index(i, j, node->M_block, node->N_block)] = node->block_dense[spamm_dense_index(i, j, node->M_block, node->N_block)];
        }
      }
    }
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
