#include "spamm.h"
#include "config.h"
#include <assert.h>
#include <math.h>
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
 * @param tolerance The accuracy target for the matrix product.
 * @param alpha The scalar factor multiplying A*B.
 * @param A_node The node of matrix A.
 * @param B_node The node of matrix B.
 * @param C_node The node of matrix C.
 * @param multiply_stream The multiply stream.
 *
 * @return The number of block matrix products.
 */
unsigned int
spamm_multiply_node (const enum spamm_multiply_algorithm_t algorithm,
    const floating_point_t tolerance,
    const floating_point_t alpha, struct spamm_node_t *A_node,
    struct spamm_node_t *B_node, struct spamm_node_t **C_node,
    struct spamm_ll_t *multiply_stream)
{
  floating_point_t beta = 1.0;
  int i, j, k, l, m;
  char bitstring_A[100];
  char bitstring_B[100];
  char bitstring_C[100];
  unsigned int mask_A, mask_B, mask_C;
  unsigned int kernel_block_N;
  struct spamm_node_t *C_child_node;
  struct spamm_multiply_stream_element_t *multiply_stream_element;
  struct spamm_ll_iterator_t *iterator_A, *iterator_B, *iterator_C;
  struct spamm_ll_node_t *linear_node_A, *linear_node_B, *linear_node_C, *linear_node_C_next;
  struct spamm_linear_quadtree_t *linear_A, *linear_B, *linear_C, *linear_C_next;
  unsigned int number_products = 0;

  /* Create new node. */
  if (*C_node == NULL)
  {
    LOG2_DEBUG("creating new C node\n");

    *C_node = spamm_new_node();

    (*C_node)->tier = A_node->tier;
    (*C_node)->tree_depth = A_node->tree_depth;

    (*C_node)->M_lower = A_node->M_lower;
    (*C_node)->M_upper = A_node->M_upper;
    (*C_node)->N_lower = B_node->N_lower;
    (*C_node)->N_upper = B_node->N_upper;

    (*C_node)->linear_tier = A_node->linear_tier;
    (*C_node)->kernel_tier = A_node->kernel_tier;
  }

  /* Recurse down the tree.
   *
   * [FIXME] This should be done in index ordering, i.e. in case we have
   * Z-curve ordering, we should recurse in that order and not simply on
   * child_{ij} with 2 nested loops.
   */
  LOG2_DEBUG("recursing...\n");
  for (i = 0; i < SPAMM_N_CHILD; ++i) {
    for (j = 0; j < SPAMM_N_CHILD; ++j) {
      for (k = 0; k < SPAMM_N_CHILD; ++k)
      {
        if (A_node->child[i][k] != NULL && B_node->child[k][j] != NULL)
        {
          if (A_node->child[i][k]->norm*B_node->child[k][j]->norm >= tolerance)
          {
            if ((*C_node)->child[i][j] == NULL)
            {
              /* Create new child node in C. */
              (*C_node)->child[i][j] = spamm_new_node();
              C_child_node = (*C_node)->child[i][j];

              C_child_node->tier = (*C_node)->tier+1;
              C_child_node->tree_depth = (*C_node)->tree_depth;

              C_child_node->M_lower = (*C_node)->M_lower+i*((*C_node)->M_upper-(*C_node)->M_lower)/SPAMM_N_CHILD;
              C_child_node->M_upper = (*C_node)->M_lower+(i+1)*((*C_node)->M_upper-(*C_node)->M_lower)/SPAMM_N_CHILD;
              C_child_node->N_lower = (*C_node)->N_lower+j*((*C_node)->N_upper-(*C_node)->N_lower)/SPAMM_N_CHILD;
              C_child_node->N_upper = (*C_node)->N_lower+(j+1)*((*C_node)->N_upper-(*C_node)->N_lower)/SPAMM_N_CHILD;

              C_child_node->linear_tier = (*C_node)->linear_tier;
              C_child_node->kernel_tier = (*C_node)->kernel_tier;

              /* Check if we are at the kernel level. */
              if (C_child_node->tier == C_child_node->kernel_tier)
              {
                /* Allocate contiguous matrix block. */
                kernel_block_N = pow(SPAMM_N_CHILD, C_child_node->tree_depth-C_child_node->kernel_tier)*SPAMM_N_BLOCK;
                C_child_node->block_dense = (floating_point_t*) spamm_allocate(sizeof(floating_point_t)*kernel_block_N*kernel_block_N);
                for (l = 0; l < kernel_block_N; ++l) {
                  for (m = 0; m < kernel_block_N; ++m)
                  {
                    C_child_node->block_dense[spamm_dense_index(l, m, kernel_block_N, kernel_block_N)] = 0.0;
                  }
                }
              }

              else if (C_child_node->tier > C_child_node->kernel_tier)
              {
                /* Point into the contiguous matrix block. */
                kernel_block_N = pow(SPAMM_N_CHILD, C_child_node->tree_depth-C_child_node->tier)*SPAMM_N_BLOCK;
                C_child_node->block_dense = (*C_node)->block_dense+kernel_block_N*kernel_block_N*(SPAMM_N_CHILD*l+k);
              }
            }

            /* Recurse. */
            number_products += spamm_multiply_node(algorithm, tolerance, alpha,
                A_node->child[i][k], B_node->child[k][j], &(*C_node)->child[i][j], multiply_stream);
          }

          else
          {
            LOG2_DEBUG("dropping product, it is below tolerance\n");
          }
        }
      }
    }
  }

  LOG2_DEBUG("done recursing, checking whether we are at the bottom.\n");
  if (A_node->tier == A_node->tree_depth && B_node->tier == B_node->tree_depth &&
      A_node->block_dense != NULL && B_node->block_dense != NULL)
  {
    /* Count this product. */
    number_products++;

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
        cublasAlloc(A_node->M_block*A_node->N_block,       sizeof(floating_point_t), &d_A);
        cublasAlloc(B_node->M_block*B_node->N_block,       sizeof(floating_point_t), &d_B);
        cublasAlloc((*C_node)->M_block*(*C_node)->N_block, sizeof(floating_point_t), &d_C);

        cublasSetMatrix(A_node->M_block,    A_node->N_block,    sizeof(floating_point_t), (void*) A_node->block_dense,    A_node->M_block,    d_A, A_node->M_block);
        cublasSetMatrix(B_node->M_block,    B_node->N_block,    sizeof(floating_point_t), (void*) B_node->block_dense,    B_node->M_block,    d_B, B_node->M_block);
        cublasSetMatrix((*C_node)->M_block, (*C_node)->N_block, sizeof(floating_point_t), (void*) (*C_node)->block_dense, (*C_node)->M_block, d_C, (*C_node)->M_block);

        cublasSgemm('N', 'N', A_node->M_block, B_node->N_block, A_node->N_block, alpha, d_A, A_node->M_block, d_B, B_node->M_block, beta, d_C, (*C_node)->M_block);
        cublasGetMatrix((*C_node)->M_block, (*C_node)->N_block, sizeof(floating_point_t), (void*) d_C, (*C_node)->M_block, (void*) (*C_node)->block_dense, (*C_node)->M_block);

        cublasFree(d_A);
        cublasFree(d_B);
        cublasFree(d_C);
#else
        spamm_sgemm_trivial('N', 'N', SPAMM_N_BLOCK, SPAMM_N_BLOCK,
            SPAMM_N_BLOCK, alpha, A_node->block_dense, SPAMM_N_BLOCK,
            B_node->block_dense, SPAMM_N_BLOCK, beta,
            (*C_node)->block_dense, SPAMM_N_BLOCK);
#endif
        break;

      case cache:
      case cache_redundant:
        /* Append this triple to the multiply stream. */
        multiply_stream_element = (struct spamm_multiply_stream_element_t*) malloc(sizeof(struct spamm_multiply_stream_element_t));
        multiply_stream_element->alpha = alpha;
        multiply_stream_element->beta = 1.0;
        multiply_stream_element->A_index = A_node->index;
        multiply_stream_element->B_index = B_node->index;
        multiply_stream_element->C_index = (*C_node)->index;
        multiply_stream_element->M_A = SPAMM_N_BLOCK;
        multiply_stream_element->N_A = SPAMM_N_BLOCK;
        multiply_stream_element->M_B = SPAMM_N_BLOCK;
        multiply_stream_element->N_B = SPAMM_N_BLOCK;
        multiply_stream_element->M_C = SPAMM_N_BLOCK;
        multiply_stream_element->N_C = SPAMM_N_BLOCK;
        multiply_stream_element->block_A_loaded_in_GPU = 0;
        multiply_stream_element->block_B_loaded_in_GPU = 0;
        multiply_stream_element->block_C_loaded_in_GPU = 0;
        multiply_stream_element->A_block_dense = A_node->block_dense;
        multiply_stream_element->B_block_dense = B_node->block_dense;
        multiply_stream_element->C_node = *C_node;
        if (algorithm == cache_redundant)
        {
          multiply_stream_element->C_block_dense = (floating_point_t*) malloc(sizeof(floating_point_t)*SPAMM_N_BLOCK*SPAMM_N_BLOCK);
        }

        else
        {
          multiply_stream_element->C_block_dense = (*C_node)->block_dense;
        }

        spamm_ll_append(multiply_stream_element, multiply_stream);
        break;

      default:
        LOG2_FATAL("unknown algorithm\n");
        exit(1);
        break;
    }
  }

  else
  {
    LOG2_DEBUG("either A or B is zero\n");
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
          linear_C = spamm_new_linear_quadtree_node(SPAMM_N_BLOCK, SPAMM_N_BLOCK, (*C_node)->linear_quadtree_memory);
          linear_C->index = mask_C;
          spamm_ll_append(linear_C, (*C_node)->linear_quadtree);

          /* Multiply blocks. */
          LOG2_DEBUG("multiplying blocks\n");
          LOG2_DEBUG("A:\n");
          if (spamm_get_loglevel() == debug) { spamm_print_dense(linear_A->M, linear_A->N, linear_A->block_dense); }
          LOG2_DEBUG("B:\n");
          if (spamm_get_loglevel() == debug) { spamm_print_dense(linear_B->M, linear_B->N, linear_B->block_dense); }

          for (i = 0; i < SPAMM_N_BLOCK; ++i) {
            for (j = 0; j < SPAMM_N_BLOCK; ++j)
            {
              linear_C->block_dense[spamm_dense_index(i, j, linear_C->M, linear_C->N)] = 0.0;

              for (k = 0; k < SPAMM_N_BLOCK; ++k)
              {
                linear_C->block_dense[spamm_dense_index(i, j, linear_C->M, linear_C->N)] +=
                  alpha*linear_A->block_dense[spamm_dense_index(i, k, SPAMM_N_BLOCK, SPAMM_N_BLOCK)]
                  *linear_B->block_dense[spamm_dense_index(k, j, SPAMM_N_BLOCK, SPAMM_N_BLOCK)];
              }
            }
          }
          LOG2_DEBUG("C:\n");
          if (spamm_get_loglevel() == debug) { spamm_print_dense(linear_C->M, linear_C->N, linear_C->block_dense); }
        }
      }
      spamm_ll_iterator_delete(&iterator_B);
    }
    spamm_ll_iterator_delete(&iterator_A);

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
          spamm_ll_delete_node(NULL, linear_node_C->next, (*C_node)->linear_quadtree);

          LOG2_DEBUG("linear_C after:\n");
          if (spamm_get_loglevel() == debug) { spamm_print_dense(linear_C->M, linear_C->N, linear_C->block_dense); }
        }

        else { break; }
      }
    }
    spamm_ll_iterator_delete(&iterator_C);

    LOG_INFO("done, C now has %u elements\n", (*C_node)->linear_quadtree->number_elements);
  }

  return number_products;
}

/** Go through the multiply stream and calculate product.
 *
 * @param cache_length Determines the number of blocks kept in the GPU.
 * @param multiply_stream The multiply stream.
 */
void
spamm_multiply_stream (const unsigned int cache_length, const struct spamm_ll_t *multiply_stream)
{
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
      if (nodedata->block_A_loaded_in_GPU == 1)
      {
        nodedata->block_A_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
        cublasFree(nodedata->A_node->device_pointer);
#endif
      }

      if (nodedata->block_B_loaded_in_GPU == 1)
      {
        nodedata->block_B_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
        cublasFree(nodedata->B_node->device_pointer);
#endif
      }

      if (nodedata->block_C_loaded_in_GPU == 1)
      {
        nodedata->block_C_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
        cublasGetMatrix(nodedata->C_node->M_block, nodedata->C_node->N_block, sizeof(floating_point_t),
            (void*) nodedata->C_node->device_pointer, nodedata->C_node->M_block,
            (void*) nodedata->C_node->block_dense, nodedata->C_node->M_block);
        cublasFree(nodedata->C_node->device_pointer);
#endif
      }

      tailnode = tailnode->next;
      head_tail_distance--;
    }

    nodedata = (struct spamm_multiply_stream_element_t*) node->data;
    if (nodedata->block_A_loaded_in_GPU == 0)
    {
      number_A_blocks_loaded++;
      nodedata->block_A_loaded_in_GPU = 1;
#ifdef HAVE_CUDA
      cublasAlloc(nodedata->A_node->M_block*nodedata->A_node->N_block, sizeof(floating_point_t), &(nodedata->A_node->device_pointer));
      cublasSetMatrix(nodedata->A_node->M_block, nodedata->A_node->N_block, sizeof(floating_point_t),
          (void*) nodedata->A_node->block_dense, nodedata->A_node->M_block,
          nodedata->A_node->device_pointer, nodedata->A_node->M_block);
#endif
    }

    if (nodedata->block_B_loaded_in_GPU == 0)
    {
      number_B_blocks_loaded++;
      nodedata->block_B_loaded_in_GPU = 1;
#ifdef HAVE_CUDA
      cublasAlloc(nodedata->B_node->M_block*nodedata->B_node->N_block, sizeof(floating_point_t), &(nodedata->B_node->device_pointer));
      cublasSetMatrix(nodedata->B_node->M_block, nodedata->B_node->N_block, sizeof(floating_point_t),
          (void*) nodedata->B_node->block_dense, nodedata->B_node->M_block,
          nodedata->B_node->device_pointer, nodedata->B_node->M_block);
#endif
    }

    if (nodedata->block_C_loaded_in_GPU == 0)
    {
      number_C_blocks_loaded++;
      nodedata->block_C_loaded_in_GPU = 1;
#ifdef HAVE_CUDA
      cublasAlloc(nodedata->C_node->M_block*nodedata->C_node->N_block, sizeof(floating_point_t), &(nodedata->C_node->device_pointer));
      cublasSetMatrix(nodedata->C_node->M_block, nodedata->C_node->N_block, sizeof(floating_point_t),
          (void*) nodedata->C_node->block_dense, nodedata->C_node->M_block,
          nodedata->C_node->device_pointer, nodedata->C_node->M_block);
#endif
    }

#if defined(HAVE_CUDA)
    cublasSgemm('N', 'N', nodedata->A_node->M_block, nodedata->B_node->N_block, nodedata->A_node->N_block,
        nodedata->alpha, nodedata->A_node->device_pointer, nodedata->A_node->M_block, nodedata->B_node->device_pointer, nodedata->B_node->M_block,
        nodedata->beta, nodedata->C_node->device_pointer, nodedata->C_node->M_block);
#elif defined(DGEMM)
    DGEMM("N", "N", &(nodedata->M_A), &(nodedata->N_B), &(nodedata->N_A),
        &(nodedata->alpha), nodedata->A_block_dense, &(nodedata->M_A),
        nodedata->B_block_dense, &(nodedata->M_B), &(nodedata->beta),
        nodedata->C_block_dense, &(nodedata->M_C));
#else
    spamm_sgemm_trivial('N', 'N', nodedata->M_A, nodedata->N_B, nodedata->N_A,
        nodedata->alpha, nodedata->A_block_dense, nodedata->M_A,
        nodedata->B_block_dense, nodedata->M_B, nodedata->beta,
        nodedata->C_block_dense, nodedata->M_C);
#endif

    head_tail_distance++;
  }

  /* Unload the remaining C bocks. */
  while (tailnode != NULL)
  {
    nodedata = (struct spamm_multiply_stream_element_t*) tailnode->data;
    if (nodedata->block_A_loaded_in_GPU == 1)
    {
      nodedata->block_A_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
      cublasFree(nodedata->A_node->device_pointer);
#endif
    }

    if (nodedata->block_B_loaded_in_GPU == 1)
    {
      nodedata->block_B_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
      cublasFree(nodedata->B_node->device_pointer);
#endif
    }

    if (nodedata->block_C_loaded_in_GPU == 1)
    {
      nodedata->block_C_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
      cublasGetMatrix(nodedata->C_node->M_block, nodedata->C_node->N_block, sizeof(floating_point_t), (void*) nodedata->C_node->device_pointer,
          nodedata->C_node->M_block, (void*) nodedata->C_node->block_dense, nodedata->C_node->M_block);
      cublasFree(nodedata->C_node->device_pointer);
#endif
    }

    tailnode = tailnode->next;
  }

  /* Print statistics. */
  LOG_INFO("blocks loaded (stream length = %u, cache = %u): A = %u, B = %u, C = %u, total = %u\n",
      multiply_stream->number_elements, cache_length,
      number_A_blocks_loaded, number_B_blocks_loaded, number_C_blocks_loaded,
      number_A_blocks_loaded+number_B_blocks_loaded+number_C_blocks_loaded);
}

/** Go through the multiply stream and resum duplicate C blocks.
 *
 * @param cache_length Determines the number of blocks kept in the GPU.
 * @param multiply_stream The multiply stream.
 */
void
spamm_resum_stream (const unsigned int cache_length, struct spamm_ll_t *multiply_stream)
{
  unsigned int i;
  struct spamm_ll_iterator_t *iterator;
  struct spamm_ll_node_t *stream_node, *next_stream_node;
  struct spamm_multiply_stream_element_t *stream_element, *next_stream_element;
  struct timeval start, stop;

  /* Sort stream on C block index. */
  LOG_INFO("resum: sorting multiply stream in C (has %u elements)\n", multiply_stream->number_elements);
  //spamm_print_multiply_stream(multiply_stream);
  gettimeofday(&start, NULL);
  spamm_ll_sort_data(spamm_compare_multiply_stream_element, spamm_swap_multiply_stream, multiply_stream);
  gettimeofday(&stop, NULL);
  LOG_INFO("resum: sorting multiply stream: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);
  //spamm_print_multiply_stream(multiply_stream);

  /* Find duplicate blocks in C and sum them. */
  LOG2_INFO("resum: adding duplicate blocks in C\n");
  iterator = spamm_ll_iterator_new(multiply_stream);
  for (stream_node = spamm_ll_iterator_first(iterator); stream_node != NULL; stream_node = spamm_ll_iterator_next(iterator))
  {
    next_stream_node = stream_node->next;
    stream_element = stream_node->data;
    //LOG_INFO("stream element: index = %u, data = %f\n", stream_element->C_index, stream_element->C_block_dense[0]);
    while (next_stream_node != NULL)
    {
      next_stream_element = stream_node->next->data;
      //LOG_INFO("next stream element: index = %u, data = %f\n", next_stream_element->C_index, next_stream_element->C_block_dense[0]);

      if (next_stream_element->C_index == stream_element->C_index)
      {
        /* Sum blocks. */
        for (i = 0; i < stream_element->M_C*stream_element->N_C; ++i)
        {
          stream_element->C_block_dense[i] += next_stream_element->C_block_dense[i];
        }
        next_stream_node = next_stream_node->next;
        spamm_ll_delete_node(NULL, stream_node->next, multiply_stream);
      }

      else { break; }
    }
  }

  /* Free memory. */
  spamm_ll_iterator_delete(&iterator);
  //spamm_print_multiply_stream(multiply_stream);
}

/** Go through the multiply stream and add C blocks to C matrix.
 *
 * @param cache_length Determines the number of blocks kept in the GPU.
 * @param multiply_stream The multiply stream.
 */
void
spamm_add_stream (const unsigned int cache_length, struct spamm_ll_t *multiply_stream)
{
  unsigned int i;
  struct spamm_ll_iterator_t *iterator;
  struct spamm_ll_node_t *stream_node, *next_stream_node;
  struct spamm_multiply_stream_element_t *stream_element, *next_stream_element;

  LOG_INFO("add: adding multiply stream to C (has %u elements)\n", multiply_stream->number_elements);
  iterator = spamm_ll_iterator_new(multiply_stream);
  for (stream_node = spamm_ll_iterator_first(iterator); stream_node != NULL; stream_node = spamm_ll_iterator_next(iterator))
  {
    stream_element = stream_node->data;
    for (i = 0; i < stream_element->M_C*stream_element->N_C; ++i)
    {
      stream_element->C_node->block_dense[i] += stream_element->C_block_dense[i];
    }
  }

  /* Free memory. */
  spamm_ll_iterator_delete(&iterator);
}

/** Computes the product
 *
 * \f$C = \alpha A \times B + \beta C\f$
 *
 * @param algorithm The algorithm to use.
 * @param tolerance The accuracy target for the matrix product.
 * @param alpha The scalar factor \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param beta The scalar factor \f$\beta\f$.
 * @param C The matrix \f$C\f$.
 *
 * @return The number of block matrix products.
 *
 * \bug Can not handle multiply of trees with different depths.
 */
unsigned int
spamm_multiply (const enum spamm_multiply_algorithm_t algorithm,
    floating_point_t tolerance,
    const floating_point_t alpha, const struct spamm_t *A,
    const struct spamm_t *B, const floating_point_t beta, struct spamm_t *C)
{
  struct timeval start, stop;
  struct timeval total_start, total_stop;
  struct spamm_ll_t *multiply_stream;
  unsigned int number_products = 0;

  double max_memory;

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

  if (tolerance < 0.0)
  {
    LOG2_INFO("tolerance is negative, using absolute value\n");
    tolerance = -tolerance;
  }

  if (beta != 1.0)
  {
    /* Multiply existing C. */
    spamm_multiply_scalar(beta, C);
  }

  if ((A->root == NULL || B->root == NULL) && C->root == NULL)
  {
    /* Nothing to be done. */
    return number_products;
  }

  switch (algorithm)
  {
    case tree:
      LOG2_INFO("using tree algorithm\n");
      number_products = spamm_multiply_node(algorithm, tolerance, alpha, A->root, B->root, &(C->root), NULL);
      break;

    case cache:
      LOG2_INFO("using cache algorithm\n");
    case cache_redundant:
      if (algorithm == cache_redundant)
      {
        LOG2_INFO("using cache (redundant) algorithm\n");
      }

      max_memory = pow(SPAMM_N_CHILD, C->tree_depth)*pow(SPAMM_N_CHILD, C->tree_depth)*pow(SPAMM_N_CHILD, A->tree_depth)
        *(sizeof(struct spamm_multiply_stream_element_t)+SPAMM_N_BLOCK*SPAMM_N_BLOCK*sizeof(floating_point_t));
      if (max_memory < 1024)
      {
        LOG_INFO("max memory usage for multiply stream: %1.2f bytes\n", max_memory);
      }

      else if (max_memory < 1024*1024)
      {
        LOG_INFO("max memory usage for multiply stream: %1.2f kB\n", max_memory/1024.);
      }

      else if (max_memory < 1024*1024*1024)
      {
        LOG_INFO("max memory usage for multiply stream: %1.2f MB\n", max_memory/1024./1024.);
      }

      else
      {
        LOG_INFO("max memory usage for multiply stream: %1.2f GB\n", max_memory/1024./1024./1024.);
      }

      gettimeofday(&total_start, NULL);
      multiply_stream = spamm_ll_new();

      gettimeofday(&start, NULL);
      number_products = spamm_multiply_node(algorithm, tolerance, alpha, A->root, B->root, &(C->root), multiply_stream);
      gettimeofday(&stop, NULL);
      LOG_INFO("symbolic multiply: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

      gettimeofday(&start, NULL);
      spamm_multiply_stream(multiply_stream->number_elements, multiply_stream);
      gettimeofday(&stop, NULL);
      LOG_INFO("stream multiply: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

      //gettimeofday(&start, NULL);
      //spamm_resum_stream(multiply_stream->number_elements, multiply_stream);
      //gettimeofday(&stop, NULL);
      //LOG_INFO("resum stream: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

      if (algorithm == cache_redundant)
      {
        gettimeofday(&start, NULL);
        spamm_add_stream(multiply_stream->number_elements, multiply_stream);
        gettimeofday(&stop, NULL);
        LOG_INFO("add stream: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);
      }

      spamm_ll_delete(spamm_delete_multiply_stream_element, &multiply_stream);

      gettimeofday(&total_stop, NULL);
      LOG_INFO("total time for multiply: %f s\n", (total_stop.tv_sec-total_start.tv_sec)+(total_stop.tv_usec-total_start.tv_usec)/(double) 1e6);
      break;

    default:
      LOG2_FATAL("unknown algorithm\n");
      exit(1);
      break;
  }

  return number_products;
}
