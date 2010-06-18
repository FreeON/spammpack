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
    LOG("reached tier %u, packing\n", node->tier);
    LOG("linear_tree_memory at %p, linear_tree at %p\n", linear_tree_memory, linear_tree);
#endif
  }

  /* Recurse more. */
  if (node->child != NULL)
  {
#ifdef SPAMM_DEBUG
    LOG("linear_tree at %p\n", linear_tree);
#endif
    for (i = 0; i < node->M_child; ++i) {
      for (j = 0; j < node->N_child; ++j)
      {
        spamm_tree_pack_subtree(linear_tier, chunksize, mask, linear_tree, linear_tree_memory, node->child[spamm_dense_index(i, j, node->M_child, node->N_child)]);
      }
    }
  }

  else
  {
#ifdef SPAMM_DEBUG
    LOG("linear_tree at %p\n", linear_tree);
#endif
    if (node->block_dense != NULL)
    {
      /* Copy data into linear tree. */
#ifdef SPAMM_DEBUG
      LOG("tier %u, storing datablock with index %u\n", node->tier, node->index);
#endif
      /* Allocate new linear quadtree entry. The memory required is the size
       * of the struct plus the size of the dense matrix block plus some
       * padding.
       */
      linear_block = (struct spamm_linear_quadtree_t*) spamm_mm_allocate(sizeof(struct spamm_linear_quadtree_t)+node->M_block*node->N_block*sizeof(float_t)+8, linear_tree_memory);
      spamm_ll_append(linear_block, linear_tree);

#ifdef SPAMM_DEBUG
      LOG("linear_block at %p\n", linear_block);
      spamm_ll_print(NULL, linear_tree);
#endif

      /* Set the linear_block->block_dense pointer to just after the pointer
       * itself. Since we allocated enough space in spamm_mm_allocate() to fit
       * the pointer and the dense matrix block, we can do that without having
       * to allocate the dense data block separately.
       */
      linear_block->block_dense = (float_t*) (((void*) linear_block)+sizeof(struct spamm_linear_quadtree_t));
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
#ifdef SPAMM_DEBUG
    LOG2("done packing, deleting original tree\n");
#endif
    if (node->block_dense != NULL)
    {
      free(node->block_dense);
      node->block_dense = NULL;
    }

    if (node->child != NULL)
    {
      for (i = 0; i < node->M_child; ++i) {
        for (j = 0; j < node->N_child; ++j)
        {
          spamm_delete_node(&node->child[spamm_dense_index(i, j, node->M_child, node->N_child)]);
        }
      }
      free(node->child);
      node->child = NULL;
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
    //LOG("tree depth is only %u, while linear_tier is %u\n", A->tree_depth, linear_tier);
    return;
  }
  A->linear_tier = linear_tier;

#ifdef SPAMM_DEBUG
  LOG2("packing tree\n");
#endif

  /* Recurse tree and convert to linear trees any subtrees below linear_tier.
   */
  if (A->root != NULL)
  {
    spamm_tree_pack_subtree(linear_tier, chunksize, mask, NULL, NULL, A->root);
  }

  else
  {
    //LOG2("A is empty, nothing to pack\n");
  }
}
