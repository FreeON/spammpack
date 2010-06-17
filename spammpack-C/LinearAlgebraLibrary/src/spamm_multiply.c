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
 * C_node = alpha*A_node*B_node + beta*C_node
 *
 * @param algorithm The algorithm to use.
 * @param alpha The scalar factor multiplying A*B.
 * @param A_node The node of matrix A.
 * @param B_node The node of matrix B.
 * @param beta The scalar factor multiplying C.
 * @param C_node The node of matrix C.
 * @param multiply_stream The multiply stream.
 */
void
spamm_multiply_node (const enum spamm_multiply_algorithm_t algorithm,
    const float_t alpha, struct spamm_node_t *A_node,
    struct spamm_node_t *B_node, const float_t beta,
    struct spamm_node_t **C_node, struct spamm_ll_t *multiply_stream)
{
  int i, j, k;
  struct spamm_node_t *C_child_node;
  struct spamm_multiply_stream_element_t *multiply_stream_element;

  if (*C_node == NULL)
  {
    /* Create new node. */
    spamm_new_node(C_node);

    (*C_node)->tier = A_node->tier+1;
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
          spamm_new_node(&((*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]));

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
            spamm_log("bug?\n", __FILE__, __LINE__);
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
          spamm_multiply_node(algorithm, alpha, A_node->child[spamm_dense_index(i, k, A_node->M_child, A_node->N_child)],
              B_node->child[spamm_dense_index(k, j, B_node->M_child, B_node->N_child)], beta,
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
            for (k = 0; k < A_node->M_block; ++k)
            {
              (*C_node)->block_dense[spamm_dense_index(i, j, (*C_node)->M_block, (*C_node)->N_block)]
                = alpha*A_node->block_dense[spamm_dense_index(i, k, A_node->M_block, A_node->N_block)]*B_node->block_dense[spamm_dense_index(k, j, B_node->M_block, B_node->N_block)]
                + beta*(*C_node)->block_dense[spamm_dense_index(i, j, (*C_node)->M_block, (*C_node)->N_block)];
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
        spamm_log("unknown algorithm\n", __FILE__, __LINE__);
        exit(1);
        break;
    }
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
  LOG("blocks loaded: stream length = %3u cache = %3u A = %3u B = %3u C = %3u total = %4u\n",
      multiply_stream->number_elements, cache_length,
      number_A_blocks_loaded, number_B_blocks_loaded, number_C_blocks_loaded,
      number_A_blocks_loaded+number_B_blocks_loaded+number_C_blocks_loaded);
}

/** Computes the product
 *
 * C = alpha*A*B + beta*C
 *
 * @param algorithm The algorithm to use.
 * @param alpha The scalar factor multiplying A*B.
 * @param A The matrix A.
 * @param B The matrix B.
 * @param beta The scalar factor multiplying C.
 * @param C The matrix C.
 *
 * \bug A pre-existing C matrix will not work right now.
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
    LOG("matrix size mismatch, A->N = %i, B->M = %i\n", A->N, B->M);
    exit(1);
  }

  if (A->M != C->M)
  {
    LOG("matrix size mismatch, A->M = %i, C->M = %i\n", A->M, C->M);
    exit(1);
  }

  if (B->N != C->N)
  {
    LOG("matrix size mismatch, B->N = %i, C->N = %i\n", B->N, C->N);
    exit(1);
  }

  if (A->N_child != B->M_child)
  {
    LOG("matrix child size mismatch, A->N_child = %i, B_child->M = %i\n", A->N_child, B->M_child);
    exit(1);
  }

  if (A->M_child != C->M_child)
  {
    LOG("matrix child size mismatch, A->M_child = %i, C->M_child = %i\n", A->M_child, C->M_child);
    exit(1);
  }

  if (B->N_child != C->N_child)
  {
    LOG("matrix child size mismatch, B->N_child = %i, C->N_child = %i\n", B->N_child, C->N_child);
    exit(1);
  }

  if (A->N_block != B->M_block)
  {
    LOG("matrix block size mismatch, A->N_block = %i, B_block->M = %i\n", A->N_block, B->M_block);
    exit(1);
  }

  if (A->M_block != C->M_block)
  {
    LOG("matrix block size mismatch, A->M_block = %i, C->M_block = %i\n", A->M_block, C->M_block);
    exit(1);
  }

  if (B->N_block != C->N_block)
  {
    LOG("matrix block size mismatch, B->N_block = %i, C->N_block = %i\n", B->N_block, C->N_block);
    exit(1);
  }

  if (C->root != NULL)
  {
    spamm_log("[FIXME] can not handle pre-existing C\n", __FILE__, __LINE__);
    exit(1);
  }

  if (beta != 1.0)
  {
    LOG("[FIXME] can not handle (beta = %e) != 1.0\n", beta);
    exit(1);
  }

  if ((A->root == NULL || B->root == NULL) && C->root == NULL)
  {
    /* Nothing to be done. */
    return;
  }

  else if ((A->root == NULL || B->root == NULL) && C->root != NULL)
  {
    spamm_add_node(0.0, NULL, beta, &(C->root));
  }

  switch (algorithm)
  {
    case tree:
      spamm_multiply_node(algorithm, alpha, A->root, B->root, beta, &(C->root), multiply_stream);
      break;

    case cache:
      multiply_stream = spamm_ll_new();

      gettimeofday(&start, NULL);
      spamm_multiply_node(algorithm, alpha, A->root, B->root, beta, &(C->root), multiply_stream);
      gettimeofday(&stop, NULL);
      LOG("tree recursion: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

      gettimeofday(&start, NULL);
      spamm_multiply_stream(30000, multiply_stream);
      gettimeofday(&stop, NULL);
      LOG("stream multiply: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

      spamm_ll_delete(multiply_stream);
      break;

    default:
      spamm_log("unknown algorithm\n", __FILE__, __LINE__);
      exit(1);
      break;
  }
}
