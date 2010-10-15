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
    struct spamm_mm_t *linear_tree_memory, struct spamm_node_t *node)
{
  int i, j;
  struct spamm_linear_quadtree_t *linear_block;

  assert(node != NULL);

  node->linear_tier = linear_tier;

  if (node->tier == linear_tier)
  {
    /* Pack. */
    node->linear_quadtree_memory = spamm_mm_new(chunksize);
    node->linear_quadtree = spamm_ll_new();
    linear_tree_memory = node->linear_quadtree_memory;
    linear_tree = node->linear_quadtree;
#ifdef SPAMM_DEBUG
    LOG_DEBUG("reached tier %u, packing\n", node->tier);
    LOG_DEBUG("linear_tree_memory at %p, linear_tree at %p\n", linear_tree_memory, linear_tree);
#endif
  }

  /* Recurse more. */
  if (node->child != NULL)
  {
#ifdef SPAMM_DEBUG
    LOG_DEBUG("linear_tree at %p\n", linear_tree);
#endif
    for (i = 0; i < SPAMM_M_CHILD; ++i) {
      for (j = 0; j < SPAMM_N_CHILD; ++j)
      {
        spamm_tree_pack_subtree(linear_tier, chunksize, mask, linear_tree, linear_tree_memory, node->child[i][j]);
      }
    }
  }

  else
  {
#ifdef SPAMM_DEBUG
    LOG_DEBUG("linear_tree at %p\n", linear_tree);
#endif
    if (node->block_dense != NULL)
    {
      /* Copy data into linear tree. */
#ifdef SPAMM_DEBUG
      LOG_DEBUG("tier %u, storing datablock with index %u\n", node->tier, node->index);
#endif
      /* Allocate new linear quadtree entry. The memory required is the size
       * of the struct plus the size of the dense matrix block plus some
       * padding.
       */
      linear_block = (struct spamm_linear_quadtree_t*) spamm_mm_allocate(sizeof(struct spamm_linear_quadtree_t)+SPAMM_M_BLOCK*SPAMM_N_BLOCK*sizeof(floating_point_t)+8, linear_tree_memory);
      spamm_ll_append(linear_block, linear_tree);

#ifdef SPAMM_DEBUG
      LOG_DEBUG("linear_block at %p\n", linear_block);
      spamm_ll_print(NULL, linear_tree);
#endif

      /* Set the linear_block->block_dense pointer to just after the pointer
       * itself. Since we allocated enough space in spamm_mm_allocate() to fit
       * the pointer and the dense matrix block, we can do that without having
       * to allocate the dense data block separately.
       */
      linear_block->block_dense = (floating_point_t*) (((void*) linear_block)+sizeof(struct spamm_linear_quadtree_t));
      linear_block->index = node->index;
      linear_block->M = SPAMM_M_BLOCK;
      linear_block->N = SPAMM_N_BLOCK;

      /* Copy data. */
      for (i = 0; i < SPAMM_M_BLOCK; ++i) {
        for (j = 0; j < SPAMM_N_BLOCK; ++j)
        {
          linear_block->block_dense[spamm_dense_index(i, j, SPAMM_M_BLOCK, SPAMM_N_BLOCK)] = node->block_dense[spamm_dense_index(i, j, SPAMM_M_BLOCK, SPAMM_N_BLOCK)];
        }
      }
    }
  }

  if (node->tier == linear_tier)
  {
    /* Delete tree, which is replaced now by linear quadtree. */
#ifdef SPAMM_DEBUG
    LOG2_DEBUG("done packing, deleting original tree\n");
#endif
    if (node->block_dense != NULL)
    {
      free(node->block_dense);
      node->block_dense = NULL;
    }

    if (node->child != NULL)
    {
      for (i = 0; i < SPAMM_M_CHILD; ++i) {
        for (j = 0; j < SPAMM_M_CHILD; ++j)
        {
          spamm_delete_node(&node->child[i][j]);
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
    LOG_DEBUG("tree depth is only %u, while linear_tier is %u\n", A->tree_depth, linear_tier);
    return;
  }
  A->linear_tier = linear_tier;

#ifdef SPAMM_DEBUG
  LOG2_DEBUG("packing tree\n");
#endif

  /* Recurse tree and convert to linear trees any subtrees below linear_tier.
   */
  if (A->root != NULL)
  {
    spamm_tree_pack_subtree(linear_tier, chunksize, mask, NULL, NULL, A->root);
  }

  else
  {
    LOG2_DEBUG("A is empty, nothing to pack\n");
  }
}
