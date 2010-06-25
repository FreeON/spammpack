#include "spamm.h"
#include "config.h"
#include <assert.h>
#include <stdlib.h>
#include <sys/time.h>

#if defined(HAVE_CUDA)
#include <cublas.h>
#endif

/** \private Computes
 *
 * C_node = alpha*A_node*B_node + C_node
 *
 * @param algorithm The algorithm to use.
 * @param alpha The scalar factor multiplying A*B.
 * @param A_node The node of matrix A.
 * @param B_node The node of matrix B.
 * @param C_node The node of matrix C.
 * @param multiply_stream The multiply stream.
 */
void
spamm_multiply_node (const enum spamm_multiply_algorithm_t algorithm,
    const float_t alpha, struct spamm_node_t *A_node,
    struct spamm_node_t *B_node, struct spamm_node_t **C_node,
    struct spamm_ll_t *multiply_stream)
{
  float_t beta = 1.0;
  int i, j, k;
  char bitstring_A[100];
  char bitstring_B[100];
  char bitstring_C[100];
  unsigned int mask_A, mask_B, mask_C;
  unsigned int bit_A, bit_B, match;
  struct spamm_node_t *C_child_node;
  struct spamm_multiply_stream_element_t *multiply_stream_element;
  struct spamm_ll_iterator_t *iterator_A, *iterator_B, *iterator_C;
  struct spamm_ll_node_t *linear_node_A, *linear_node_B, *linear_node_C, *linear_node_C_next;
  struct spamm_linear_quadtree_t *linear_A, *linear_B, *linear_C, *linear_C_next;

  /* Create new node. */
  if (*C_node == NULL)
  {
    *C_node = spamm_new_node();

    (*C_node)->tier = A_node->tier;
    (*C_node)->tree_depth = A_node->tree_depth;

    (*C_node)->M_lower = A_node->M_lower;
    (*C_node)->M_upper = A_node->M_upper;
    (*C_node)->N_lower = B_node->N_lower;
    (*C_node)->N_upper = B_node->N_upper;

    (*C_node)->M_child = A_node->M_child;
    (*C_node)->N_child = B_node->N_child;

    (*C_node)->threshold = A_node->threshold;

    (*C_node)->linear_tier = A_node->linear_tier;

    (*C_node)->M_block = A_node->M_block;
    (*C_node)->N_block = B_node->N_block;
  }

  if (A_node->child != NULL && B_node->child != NULL)
  {
    if ((*C_node)->child == NULL)
    {
      /* Allocate children nodes. */
      (*C_node)->child = (struct spamm_node_t**) malloc(sizeof(struct spamm_node_t*)*(*C_node)->M_child*(*C_node)->N_child);

      /* Create children nodes. */
      for (i = 0; i < (*C_node)->M_child; ++i) {
        for (j = 0; j < (*C_node)->N_child; ++j)
        {
          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)] = spamm_new_node();

          C_child_node = (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)];
          C_child_node->tier = (*C_node)->tier+1;
          C_child_node->tree_depth = (*C_node)->tree_depth;

          C_child_node->M_lower = (*C_node)->M_lower+i*((*C_node)->M_upper-(*C_node)->M_lower)/(*C_node)->M_child;
          C_child_node->M_upper = (*C_node)->M_lower+(i+1)*((*C_node)->M_upper-(*C_node)->M_lower)/(*C_node)->M_child;
          C_child_node->N_lower = (*C_node)->N_lower+j*((*C_node)->N_upper-(*C_node)->N_lower)/(*C_node)->N_child;
          C_child_node->N_upper = (*C_node)->N_lower+(j+1)*((*C_node)->N_upper-(*C_node)->N_lower)/(*C_node)->N_child;

          C_child_node->M_child = (*C_node)->M_child;
          C_child_node->N_child = (*C_node)->N_child;

          C_child_node->threshold = (*C_node)->threshold;

          C_child_node->linear_tier = (*C_node)->linear_tier;

          C_child_node->M_block = (*C_node)->M_block;
          C_child_node->N_block = (*C_node)->N_block;
        }
      }

      /* Space-Filling curve ordering. */
      if ((*C_node)->M_child == 2 && (*C_node)->N_child == 2)
      {
        /* Z-curve curve ordering. */
        switch ((*C_node)->ordering)
        {
          case none:
            (*C_node)->child[spamm_dense_index(0, 0, (*C_node)->M_child, (*C_node)->N_child)]->ordering = P;
            (*C_node)->child[spamm_dense_index(0, 0, (*C_node)->M_child, (*C_node)->N_child)]->index = 0;
            (*C_node)->child[spamm_dense_index(0, 1, (*C_node)->M_child, (*C_node)->N_child)]->ordering = P;
            (*C_node)->child[spamm_dense_index(0, 1, (*C_node)->M_child, (*C_node)->N_child)]->index = 1;
            (*C_node)->child[spamm_dense_index(1, 0, (*C_node)->M_child, (*C_node)->N_child)]->ordering = P;
            (*C_node)->child[spamm_dense_index(1, 0, (*C_node)->M_child, (*C_node)->N_child)]->index = 2;
            (*C_node)->child[spamm_dense_index(1, 1, (*C_node)->M_child, (*C_node)->N_child)]->ordering = P;
            (*C_node)->child[spamm_dense_index(1, 1, (*C_node)->M_child, (*C_node)->N_child)]->index = 3;
            break;

          case P:
            (*C_node)->child[spamm_dense_index(0, 0, (*C_node)->M_child, (*C_node)->N_child)]->ordering = P;
            (*C_node)->child[spamm_dense_index(0, 0, (*C_node)->M_child, (*C_node)->N_child)]->index = (*C_node)->index*4+0;
            (*C_node)->child[spamm_dense_index(0, 1, (*C_node)->M_child, (*C_node)->N_child)]->ordering = P;
            (*C_node)->child[spamm_dense_index(0, 1, (*C_node)->M_child, (*C_node)->N_child)]->index = (*C_node)->index*4+1;
            (*C_node)->child[spamm_dense_index(1, 0, (*C_node)->M_child, (*C_node)->N_child)]->ordering = P;
            (*C_node)->child[spamm_dense_index(1, 0, (*C_node)->M_child, (*C_node)->N_child)]->index = (*C_node)->index*4+2;
            (*C_node)->child[spamm_dense_index(1, 1, (*C_node)->M_child, (*C_node)->N_child)]->ordering = P;
            (*C_node)->child[spamm_dense_index(1, 1, (*C_node)->M_child, (*C_node)->N_child)]->index = (*C_node)->index*4+3;
            break;

          default:
            LOG2_FATAL("bug?\n");
            exit(1);
            break;
        }
      }
    }

    /* Recurse down the tree.
     *
     * [FIXME] This should be done in index ordering, i.e. in case we have
     * Z-curve ordering, we should recurse in that order and not simply on
     * child_{ij} with 2 nested loops.
     */
    for (i = 0; i < (*C_node)->M_child; ++i) {
      for (j = 0; j < (*C_node)->N_child; ++j) {
        for (k = 0; k < A_node->N_child; ++k)
        {
          spamm_multiply_node(algorithm, alpha,
              A_node->child[spamm_dense_index(i, k, A_node->M_child, A_node->N_child)],
              B_node->child[spamm_dense_index(k, j, B_node->M_child, B_node->N_child)],
              &(*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)],
              multiply_stream);
        }
      }
    }
  }

  if (A_node->block_dense != NULL && B_node->block_dense != NULL)
  {
    if ((*C_node)->block_dense == NULL)
    {
      /* Create empty dense block. */
      (*C_node)->block_dense = (float_t*) malloc(sizeof(float_t)*(*C_node)->M_block*(*C_node)->N_block);
      for (i = 0; i < (*C_node)->M_block; ++i) {
        for (j = 0; j < (*C_node)->N_block; ++j)
        {
          (*C_node)->block_dense[spamm_dense_index(i, j, (*C_node)->M_block, (*C_node)->N_block)] = 0;
        }
      }
    }

    switch (algorithm)
    {
      case tree:
        LOG_DEBUG("multiplying C(%u:%u,%u:%u) += %f*A(%u:%u,%u:%u)*B(%u:%u,%u:%u)\n",
            (*C_node)->M_lower, (*C_node)->M_upper-1, (*C_node)->N_lower, (*C_node)->N_upper-1,
            alpha,
            A_node->M_lower, A_node->M_upper-1, A_node->N_lower, A_node->N_upper-1,
            B_node->M_lower, B_node->M_upper-1, B_node->N_lower, B_node->N_upper-1);
#ifdef DGEMM
        DGEMM("N", "N", &(A_node->M_block), &(B_node->N_block), &(A_node->N_block),
            &alpha, A_node->block_dense, &(A_node->M_block), B_node->block_dense, &(B_node->M_block),
            &beta, (*C_node)->block_dense, &((*C_node)->M_block));
#elif defined(HAVE_CUDA)
        cublasAlloc(A_node->M_block*A_node->N_block,       sizeof(float_t), &d_A);
        cublasAlloc(B_node->M_block*B_node->N_block,       sizeof(float_t), &d_B);
        cublasAlloc((*C_node)->M_block*(*C_node)->N_block, sizeof(float_t), &d_C);

        cublasSetMatrix(A_node->M_block,    A_node->N_block,    sizeof(float_t), (void*) A_node->block_dense,    A_node->M_block,    d_A, A_node->M_block);
        cublasSetMatrix(B_node->M_block,    B_node->N_block,    sizeof(float_t), (void*) B_node->block_dense,    B_node->M_block,    d_B, B_node->M_block);
        cublasSetMatrix((*C_node)->M_block, (*C_node)->N_block, sizeof(float_t), (void*) (*C_node)->block_dense, (*C_node)->M_block, d_C, (*C_node)->M_block);

        cublasSgemm('N', 'N', A_node->M_block, B_node->N_block, A_node->N_block, alpha, d_A, A_node->M_block, d_B, B_node->M_block, beta, d_C, (*C_node)->M_block);
        cublasGetMatrix((*C_node)->M_block, (*C_node)->N_block, sizeof(float_t), (void*) d_C, (*C_node)->M_block, (void*) (*C_node)->block_dense, (*C_node)->M_block);

        cublasFree(d_A);
        cublasFree(d_B);
        cublasFree(d_C);
#else
        for (i = 0; i < (*C_node)->M_block; ++i) {
          for (j = 0; j < (*C_node)->N_block; ++j) {
            for (k = 0; k < A_node->N_block; ++k)
            {
              (*C_node)->block_dense[spamm_dense_index(i, j, (*C_node)->M_block, (*C_node)->N_block)]
                += alpha*A_node->block_dense[spamm_dense_index(i, k, A_node->M_block, A_node->N_block)]*B_node->block_dense[spamm_dense_index(k, j, B_node->M_block, B_node->N_block)];
            }
          }
        }
#endif
        break;

      case cache:
        /* Append this triple to the multiply stream. */
        multiply_stream_element = (struct spamm_multiply_stream_element_t*) malloc(sizeof(struct spamm_multiply_stream_element_t));
        multiply_stream_element->alpha = alpha;
        multiply_stream_element->beta = beta;
        multiply_stream_element->A_index = A_node->index;
        multiply_stream_element->B_index = B_node->index;
        multiply_stream_element->C_index = (*C_node)->index;
        multiply_stream_element->A_node = A_node;
        multiply_stream_element->B_node = B_node;
        multiply_stream_element->C_node = *C_node;

        spamm_ll_append(multiply_stream_element, multiply_stream);
        break;

      default:
        LOG2_FATAL("unknown algorithm\n");
        exit(1);
        break;
    }
  }

  if (A_node->linear_quadtree != NULL && B_node->linear_quadtree != NULL)
  {
    LOG2_DEBUG("entering linear tree...\n");

    /* Create new linear tree in C. */
    if ((*C_node)->linear_quadtree == NULL)
    {
      LOG2_DEBUG("creating new linear quadtree in C\n");
      (*C_node)->linear_quadtree = spamm_ll_new();
      (*C_node)->linear_quadtree_memory = spamm_mm_new(A_node->linear_quadtree_memory->chunksize);
    }

    /* Create convolution of A and B. */
    LOG2_INFO("creating convolution\n");
    iterator_A = spamm_ll_iterator_new(A_node->linear_quadtree);
    for (linear_node_A = spamm_ll_iterator_first(iterator_A); linear_node_A != NULL; linear_node_A = spamm_ll_iterator_next(iterator_A))
    {
      linear_A = linear_node_A->data;

      spamm_int_to_binary(linear_A->index, A_node->tree_depth*2, bitstring_A);
      LOG_DEBUG("looking at A: %s\n", bitstring_A);

      /* Loop through B to find matching blocks. */
      iterator_B = spamm_ll_iterator_new(B_node->linear_quadtree);
      for (linear_node_B = spamm_ll_iterator_first(iterator_B); linear_node_B != NULL; linear_node_B = spamm_ll_iterator_next(iterator_B))
      {
        linear_B = linear_node_B->data;

        spamm_int_to_binary(linear_A->index, A_node->tree_depth*2, bitstring_A);
        spamm_int_to_binary(linear_B->index, B_node->tree_depth*2, bitstring_B);

        LOG_DEBUG("trying to match A: %s and B: %s\n", bitstring_A, bitstring_B);

        /* We only need to compare the tiers below linear_tier. */
        mask_A = spamm_mask(linear_A->index, A_node->tree_depth*2, j_mask);
        mask_B = spamm_mask(linear_B->index, B_node->tree_depth*2, i_mask);

        spamm_int_to_binary(mask_A, A_node->tree_depth*2, bitstring_A);
        spamm_int_to_binary(mask_B, B_node->tree_depth*2, bitstring_B);

        LOG_DEBUG("mask_j(A) = %s and mask_i(B) = %s\n", bitstring_A, bitstring_B);

        if (mask_A == (mask_B >> 1))
        {
          mask_A = spamm_mask(linear_A->index, A_node->tree_depth*2, i_mask);
          mask_B = spamm_mask(linear_B->index, B_node->tree_depth*2, j_mask);
          mask_C = mask_A | mask_B;
          spamm_int_to_binary(mask_A, A_node->tree_depth*2, bitstring_A);
          spamm_int_to_binary(mask_B, B_node->tree_depth*2, bitstring_B);
          spamm_int_to_binary(mask_C, A_node->tree_depth*2, bitstring_C);

          LOG_DEBUG("found matching blocks in A and B: mask_i(A) = %s, mask_j(B) = %s, C = %s\n", bitstring_A, bitstring_B, bitstring_C);

          /* Create new block for C. */
          linear_C = spamm_new_linear_quadtree_node((*C_node)->M_block, (*C_node)->N_block, (*C_node)->linear_quadtree_memory);
          linear_C->index = mask_C;
          spamm_ll_append(linear_C, (*C_node)->linear_quadtree);

          /* Multiply blocks. */
          LOG2_DEBUG("multiplying blocks\n");
          LOG2_DEBUG("A:\n");
          if (spamm_get_loglevel() == debug) { spamm_print_dense(linear_A->M, linear_A->N, linear_A->block_dense); }
          LOG2_DEBUG("B:\n");
          if (spamm_get_loglevel() == debug) { spamm_print_dense(linear_B->M, linear_B->N, linear_B->block_dense); }

          for (i = 0; i < A_node->M_block; ++i) {
            for (j = 0; j < B_node->N_block; ++j)
            {
              linear_C->block_dense[spamm_dense_index(i, j, (*C_node)->M_block, (*C_node)->N_block)] = 0.0;

              for (k = 0; k < A_node->N_block; ++k)
              {
                linear_C->block_dense[spamm_dense_index(i, j, (*C_node)->M_block, (*C_node)->N_block)] +=
                  alpha*linear_A->block_dense[spamm_dense_index(i, k, A_node->M_block, A_node->N_block)]
                  *linear_B->block_dense[spamm_dense_index(k, j, B_node->M_block, B_node->N_block)];
              }
            }
          }
          LOG2_DEBUG("C:\n");
          if (spamm_get_loglevel() == debug) { spamm_print_dense(linear_C->M, linear_C->N, linear_C->block_dense); }
        }
      }
    }

    /* Sort blocks in C. */
    LOG_INFO("sorting linear tree in C (has %u elements)\n", (*C_node)->linear_quadtree->number_elements);
    spamm_ll_sort_data(spamm_compare_int, spamm_swap_linear_quadtree, (*C_node)->linear_quadtree);

    /* Find duplicate blocks in C and sum them. */
    LOG2_INFO("adding duplicate blocks in C\n");
    iterator_C = spamm_ll_iterator_new((*C_node)->linear_quadtree);
    for (linear_node_C = spamm_ll_iterator_first(iterator_C); linear_node_C != NULL; linear_node_C = spamm_ll_iterator_next(iterator_C))
    {
      linear_node_C_next = linear_node_C->next;

      linear_C = linear_node_C->data;
      while (linear_node_C_next != NULL)
      {
        linear_C_next = linear_node_C->next->data;

        if (linear_C_next->index == linear_C->index)
        {
          /* Sum blocks. */
          spamm_int_to_binary(linear_C->index, (*C_node)->tree_depth*2, bitstring_C);
          LOG_DEBUG("found 2 C blocks to sum: index = %s\n", bitstring_C);
          LOG2_DEBUG("linear_C before:\n");
          if (spamm_get_loglevel() == debug) { spamm_print_dense(linear_C->M, linear_C->N, linear_C->block_dense); }
          LOG2_DEBUG("linear_C_next:\n");
          if (spamm_get_loglevel() == debug) { spamm_print_dense(linear_C_next->M, linear_C_next->N, linear_C_next->block_dense); }
          for (i = 0; i < linear_C->M*linear_C->N; ++i)
          {
            linear_C->block_dense[i] += linear_C_next->block_dense[i];
          }
          linear_node_C_next = linear_node_C_next->next;
          spamm_ll_delete_node(linear_node_C->next, (*C_node)->linear_quadtree);

          LOG2_DEBUG("linear_C after:\n");
          if (spamm_get_loglevel() == debug) { spamm_print_dense(linear_C->M, linear_C->N, linear_C->block_dense); }
        }

        else { break; }
      }
    }

    LOG_INFO("done, C now has %u elements\n", (*C_node)->linear_quadtree->number_elements);
  }
}

/** Go through the multiply stream and calculate product.
 *
 * @param cache_length Determines the number of blocks kept in the GPU.
 * @param multiply_stream The multiply stream.
 */
void
spamm_multiply_stream (const unsigned int cache_length, const struct spamm_ll_t *multiply_stream)
{
#if ! defined(HAVE_CUDA) && ! defined(DGEMM)
  int i, j, k;
#endif

  unsigned int head_tail_distance;

  unsigned int number_A_blocks_loaded;
  unsigned int number_B_blocks_loaded;
  unsigned int number_C_blocks_loaded;

  struct spamm_ll_node_t *node, *tailnode;
  struct spamm_multiply_stream_element_t *nodedata;

  assert(multiply_stream != NULL);

  number_A_blocks_loaded = 0;
  number_B_blocks_loaded = 0;
  number_C_blocks_loaded = 0;
  head_tail_distance = 0;
  tailnode = multiply_stream->first;
  for (node = multiply_stream->first; node != NULL; node = node->next)
  {
    if (head_tail_distance >= cache_length)
    {
      nodedata = (struct spamm_multiply_stream_element_t*) tailnode->data;
      if (nodedata->A_node->block_loaded_in_GPU == 1)
      {
        nodedata->A_node->block_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
        cublasFree(nodedata->A_node->device_pointer);
#endif
      }

      if (nodedata->B_node->block_loaded_in_GPU == 1)
      {
        nodedata->B_node->block_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
        cublasFree(nodedata->B_node->device_pointer);
#endif
      }

      if (nodedata->C_node->block_loaded_in_GPU == 1)
      {
        nodedata->C_node->block_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
        cublasGetMatrix(nodedata->C_node->M_block, nodedata->C_node->N_block, sizeof(float_t),
            (void*) nodedata->C_node->device_pointer, nodedata->C_node->M_block,
            (void*) nodedata->C_node->block_dense, nodedata->C_node->M_block);
        cublasFree(nodedata->C_node->device_pointer);
#endif
      }

      tailnode = tailnode->next;
      head_tail_distance--;
    }

    nodedata = (struct spamm_multiply_stream_element_t*) node->data;
    if (nodedata->A_node->block_loaded_in_GPU == 0)
    {
      number_A_blocks_loaded++;
      nodedata->A_node->block_loaded_in_GPU = 1;
#ifdef HAVE_CUDA
      cublasAlloc(nodedata->A_node->M_block*nodedata->A_node->N_block, sizeof(float_t), &(nodedata->A_node->device_pointer));
      cublasSetMatrix(nodedata->A_node->M_block, nodedata->A_node->N_block, sizeof(float_t),
          (void*) nodedata->A_node->block_dense, nodedata->A_node->M_block,
          nodedata->A_node->device_pointer, nodedata->A_node->M_block);
#endif
    }

    if (nodedata->B_node->block_loaded_in_GPU == 0)
    {
      number_B_blocks_loaded++;
      nodedata->B_node->block_loaded_in_GPU = 1;
#ifdef HAVE_CUDA
      cublasAlloc(nodedata->B_node->M_block*nodedata->B_node->N_block, sizeof(float_t), &(nodedata->B_node->device_pointer));
      cublasSetMatrix(nodedata->B_node->M_block, nodedata->B_node->N_block, sizeof(float_t),
          (void*) nodedata->B_node->block_dense, nodedata->B_node->M_block,
          nodedata->B_node->device_pointer, nodedata->B_node->M_block);
#endif
    }

    if (nodedata->C_node->block_loaded_in_GPU == 0)
    {
      number_C_blocks_loaded++;
      nodedata->C_node->block_loaded_in_GPU = 1;
#ifdef HAVE_CUDA
      cublasAlloc(nodedata->C_node->M_block*nodedata->C_node->N_block, sizeof(float_t), &(nodedata->C_node->device_pointer));
      cublasSetMatrix(nodedata->C_node->M_block, nodedata->C_node->N_block, sizeof(float_t),
          (void*) nodedata->C_node->block_dense, nodedata->C_node->M_block,
          nodedata->C_node->device_pointer, nodedata->C_node->M_block);
#endif
    }

#if defined(HAVE_CUDA)
    cublasSgemm('N', 'N', nodedata->A_node->M_block, nodedata->B_node->N_block, nodedata->A_node->N_block,
        nodedata->alpha, nodedata->A_node->device_pointer, nodedata->A_node->M_block, nodedata->B_node->device_pointer, nodedata->B_node->M_block,
        nodedata->beta, nodedata->C_node->device_pointer, nodedata->C_node->M_block);
#elif defined(DGEMM)
    DGEMM("N", "N", &(nodedata->A_node->M_block), &(nodedata->B_node->N_block), &(nodedata->A_node->N_block),
        &(nodedata->alpha), nodedata->A_node->block_dense, &(nodedata->A_node->M_block), nodedata->B_node->block_dense, &(nodedata->B_node->M_block),
        &(nodedata->beta), nodedata->C_node->block_dense, &(nodedata->C_node->M_block));
#else
    for (i = 0; i < nodedata->C_node->M_block; ++i) {
      for (j = 0; j < nodedata->C_node->N_block; ++j) {
        for (k = 0; k < nodedata->A_node->M_block; ++k)
        {
          nodedata->C_node->block_dense[spamm_dense_index(i, j, nodedata->C_node->M_block, nodedata->C_node->N_block)]
            = nodedata->alpha*nodedata->A_node->block_dense[spamm_dense_index(i, k, nodedata->A_node->M_block, nodedata->A_node->N_block)]
            *nodedata->B_node->block_dense[spamm_dense_index(k, j, nodedata->B_node->M_block, nodedata->B_node->N_block)]
            + nodedata->beta*nodedata->C_node->block_dense[spamm_dense_index(i, j, nodedata->C_node->M_block, nodedata->C_node->N_block)];
        }
      }
    }
#endif

    head_tail_distance++;
  }

  /* Unload the remaining C bocks. */
  while (tailnode != NULL)
  {
    nodedata = (struct spamm_multiply_stream_element_t*) tailnode->data;
    if (nodedata->A_node->block_loaded_in_GPU == 1)
    {
      nodedata->A_node->block_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
      cublasFree(nodedata->A_node->device_pointer);
#endif
    }

    if (nodedata->B_node->block_loaded_in_GPU == 1)
    {
      nodedata->B_node->block_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
      cublasFree(nodedata->B_node->device_pointer);
#endif
    }

    if (nodedata->C_node->block_loaded_in_GPU == 1)
    {
      nodedata->C_node->block_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
      cublasGetMatrix(nodedata->C_node->M_block, nodedata->C_node->N_block, sizeof(float_t), (void*) nodedata->C_node->device_pointer,
          nodedata->C_node->M_block, (void*) nodedata->C_node->block_dense, nodedata->C_node->M_block);
      cublasFree(nodedata->C_node->device_pointer);
#endif
    }

    tailnode = tailnode->next;
  }

  /* Print statistics. */
  LOG_DEBUG("blocks loaded: stream length = %3u cache = %3u A = %3u B = %3u C = %3u total = %4u\n",
      multiply_stream->number_elements, cache_length,
      number_A_blocks_loaded, number_B_blocks_loaded, number_C_blocks_loaded,
      number_A_blocks_loaded+number_B_blocks_loaded+number_C_blocks_loaded);
}

/** Computes the product
 *
 * \f$C = \alpha A \times B + \beta C\f$
 *
 * @param algorithm The algorithm to use.
 * @param alpha The scalar factor \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param beta The scalar factor \f$\beta\f$.
 * @param C The matrix \f$C\f$.
 *
 * \bug Can not handle multiply of trees with different depths.
 */
void
spamm_multiply (const enum spamm_multiply_algorithm_t algorithm,
    const float_t alpha, const struct spamm_t *A,
    const struct spamm_t *B, const float_t beta, struct spamm_t *C)
{
  struct timeval start, stop;
  struct spamm_ll_t *multiply_stream;

  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

  if (A->N != B->M)
  {
    LOG_FATAL("matrix size mismatch, A->N = %i, B->M = %i\n", A->N, B->M);
    exit(1);
  }

  if (A->M != C->M)
  {
    LOG_FATAL("matrix size mismatch, A->M = %i, C->M = %i\n", A->M, C->M);
    exit(1);
  }

  if (B->N != C->N)
  {
    LOG_FATAL("matrix size mismatch, B->N = %i, C->N = %i\n", B->N, C->N);
    exit(1);
  }

  if (A->N_child != B->M_child)
  {
    LOG_FATAL("matrix child size mismatch, A->N_child = %i, B_child->M = %i\n", A->N_child, B->M_child);
    exit(1);
  }

  if (A->M_child != C->M_child)
  {
    LOG_FATAL("matrix child size mismatch, A->M_child = %i, C->M_child = %i\n", A->M_child, C->M_child);
    exit(1);
  }

  if (B->N_child != C->N_child)
  {
    LOG_FATAL("matrix child size mismatch, B->N_child = %i, C->N_child = %i\n", B->N_child, C->N_child);
    exit(1);
  }

  if (A->N_block != B->M_block)
  {
    LOG_FATAL("matrix block size mismatch, A->N_block = %i, B_block->M = %i\n", A->N_block, B->M_block);
    exit(1);
  }

  if (A->M_block != C->M_block)
  {
    LOG_FATAL("matrix block size mismatch, A->M_block = %i, C->M_block = %i\n", A->M_block, C->M_block);
    exit(1);
  }

  if (B->N_block != C->N_block)
  {
    LOG_FATAL("matrix block size mismatch, B->N_block = %i, C->N_block = %i\n", B->N_block, C->N_block);
    exit(1);
  }

  if (A->tree_depth != B->tree_depth || A->tree_depth != C->tree_depth)
  {
    /* The basic problem is that different tree depths means different amount
     * of padding and it is not clear to me how to recursively match up the
     * right matrix elements when the padded matrix sizes do not match. If we
     * really need this feature, a way around this would be to repad the
     * more shallow trees so that all have the same depth.
     */
    LOG2_FATAL("trees have different depths.\n");
    exit(1);
  }

  if (beta != 1.0)
  {
    /* Multiply existing C. */
    spamm_multiply_scalar(beta, C);
  }

  if ((A->root == NULL || B->root == NULL) && C->root == NULL)
  {
    /* Nothing to be done. */
    return;
  }

  switch (algorithm)
  {
    case tree:
      spamm_multiply_node(algorithm, alpha, A->root, B->root, &(C->root), multiply_stream);
      break;

    case cache:
      multiply_stream = spamm_ll_new();

      gettimeofday(&start, NULL);
      spamm_multiply_node(algorithm, alpha, A->root, B->root, &(C->root), multiply_stream);
      gettimeofday(&stop, NULL);
      LOG_DEBUG("tree recursion: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

      gettimeofday(&start, NULL);
      spamm_multiply_stream(30000, multiply_stream);
      gettimeofday(&stop, NULL);
      LOG_DEBUG("stream multiply: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

      spamm_ll_delete(free, &multiply_stream);
      break;

    default:
      LOG2_FATAL("unknown algorithm\n");
      exit(1);
      break;
  }
}
