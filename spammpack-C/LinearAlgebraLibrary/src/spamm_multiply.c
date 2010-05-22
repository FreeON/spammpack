/** @file */

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
 */
void
spamm_multiply_node (const float_t alpha, struct spamm_node_t *A_node,
    struct spamm_node_t *B_node, const float_t beta,
    struct spamm_node_t **C_node, struct spamm_multiply_stream_t *multiply_stream)
{
  int i, j, k;

  if (*C_node == NULL)
  {
    /* Create new node. */
    spamm_new_node(C_node);

    (*C_node)->tier = A_node->tier+1;

    (*C_node)->M_lower = A_node->M_lower;
    (*C_node)->M_upper = A_node->M_upper;
    (*C_node)->N_lower = B_node->N_lower;
    (*C_node)->N_upper = B_node->N_upper;

    (*C_node)->M_child = A_node->M_child;
    (*C_node)->N_child = B_node->N_child;

    (*C_node)->threshold = A_node->threshold;

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

          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->tier = (*C_node)->tier+1;

          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->M_lower = (*C_node)->M_lower+i*((*C_node)->M_upper-(*C_node)->M_lower)/(*C_node)->M_child;
          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->M_upper = (*C_node)->M_lower+(i+1)*((*C_node)->M_upper-(*C_node)->M_lower)/(*C_node)->M_child;
          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->N_lower = (*C_node)->N_lower+j*((*C_node)->N_upper-(*C_node)->N_lower)/(*C_node)->N_child;
          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->N_upper = (*C_node)->N_lower+(j+1)*((*C_node)->N_upper-(*C_node)->N_lower)/(*C_node)->N_child;

          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->M_child = (*C_node)->M_child;
          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->N_child = (*C_node)->N_child;

          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->threshold = (*C_node)->threshold;

          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->M_block = (*C_node)->M_block;
          (*C_node)->child[spamm_dense_index(i, j, (*C_node)->M_child, (*C_node)->N_child)]->N_block = (*C_node)->N_block;
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

    for (i = 0; i < (*C_node)->M_child; ++i) {
      for (j = 0; j < (*C_node)->N_child; ++j) {
        for (k = 0; k < A_node->N_child; ++k)
        {
          spamm_multiply_node(alpha, A_node->child[spamm_dense_index(i, k, A_node->M_child, A_node->N_child)],
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

    /* Append this triple to the multiply stream. */
    spamm_ll_append(alpha, beta, A_node->index, A_node, B_node->index, B_node, (*C_node)->index, *C_node, multiply_stream);
  }
}

/** Go through the multiply stream and calculate product.
 *
 * @param cache_length Determines the number of blocks kept in the GPU.
 * @param multiply_stream The multiply stream.
 */
void
spamm_multiply_stream (const unsigned int cache_length, const struct spamm_multiply_stream_t *multiply_stream)
{
  int i, j, k;

  unsigned int head_tail_distance;

  unsigned int number_A_blocks_loaded;
  unsigned int number_B_blocks_loaded;
  unsigned int number_C_blocks_loaded;

  struct spamm_multiply_stream_node_t *node, *tailnode;

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
      if (tailnode->A_node->block_loaded_in_GPU == 1)
      {
        tailnode->A_node->block_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
        cublasFree(tailnode->A_node->device_pointer);
#endif
      }

      if (tailnode->B_node->block_loaded_in_GPU == 1)
      {
        tailnode->B_node->block_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
        cublasFree(tailnode->B_node->device_pointer);
#endif
      }

      if (tailnode->C_node->block_loaded_in_GPU == 1)
      {
        tailnode->C_node->block_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
        cublasGetMatrix(tailnode->C_node->M_block, tailnode->C_node->N_block, sizeof(float_t), (void*) tailnode->C_node->device_pointer,
            tailnode->C_node->M_block, (void*) tailnode->C_node->block_dense, tailnode->C_node->M_block);
        cublasFree(tailnode->C_node->device_pointer);
#endif
      }

      tailnode = tailnode->next;
      head_tail_distance--;
    }

    if (node->A_node->block_loaded_in_GPU == 0)
    {
      number_A_blocks_loaded++;
      node->A_node->block_loaded_in_GPU = 1;
#ifdef HAVE_CUDA
      cublasAlloc(node->A_node->M_block*node->A_node->N_block, sizeof(float_t), &(node->A_node->device_pointer));
      cublasSetMatrix(node->A_node->M_block, node->A_node->N_block, sizeof(float_t), (void*) node->A_node->block_dense, node->A_node->M_block, node->A_node->device_pointer, node->A_node->M_block);
#endif
    }

    if (node->B_node->block_loaded_in_GPU == 0)
    {
      number_B_blocks_loaded++;
      node->B_node->block_loaded_in_GPU = 1;
#ifdef HAVE_CUDA
      cublasAlloc(node->B_node->M_block*node->B_node->N_block, sizeof(float_t), &(node->B_node->device_pointer));
      cublasSetMatrix(node->B_node->M_block, node->B_node->N_block, sizeof(float_t), (void*) node->B_node->block_dense, node->B_node->M_block, node->B_node->device_pointer, node->B_node->M_block);
#endif
    }

    if (node->C_node->block_loaded_in_GPU == 0)
    {
      number_C_blocks_loaded++;
      node->C_node->block_loaded_in_GPU = 1;
#ifdef HAVE_CUDA
      cublasAlloc(node->C_node->M_block*node->C_node->N_block, sizeof(float_t), &(node->C_node->device_pointer));
      cublasSetMatrix(node->C_node->M_block, node->C_node->N_block, sizeof(float_t), (void*) node->C_node->block_dense, node->C_node->M_block, node->C_node->device_pointer, node->C_node->M_block);
#endif
    }

#if defined(HAVE_CUDA)
    cublasSgemm('N', 'N', node->A_node->M_block, node->B_node->N_block, node->A_node->N_block,
        node->alpha, node->A_node->device_pointer, node->A_node->M_block, node->B_node->device_pointer, node->B_node->M_block,
        node->beta, node->C_node->device_pointer, node->C_node->M_block);
#elif defined(DGEMM)
    DGEMM("N", "N", &(node->A_node->M_block), &(node->B_node->N_block), &(node->A_node->N_block),
        &(node->alpha), node->A_node->block_dense, &(node->A_node->M_block), node->B_node->block_dense, &(node->B_node->M_block),
        &(node->beta), node->C_node->block_dense, &(node->C_node->M_block));
#else
    for (i = 0; i < node->C_node->M_block; ++i) {
      for (j = 0; j < node->C_node->N_block; ++j) {
        for (k = 0; k < node->A_node->M_block; ++k)
        {
          node->C_node->block_dense[spamm_dense_index(i, j, node->C_node->M_block, node->C_node->N_block)]
            = node->alpha*node->A_node->block_dense[spamm_dense_index(i, k, node->A_node->M_block, node->A_node->N_block)]
            *node->B_node->block_dense[spamm_dense_index(k, j, node->B_node->M_block, node->B_node->N_block)]
            + node->beta*node->C_node->block_dense[spamm_dense_index(i, j, node->C_node->M_block, node->C_node->N_block)];
        }
      }
    }
#endif

    head_tail_distance++;
  }

  /* Unload the remaining C bocks. */
  while (tailnode != NULL)
  {
    if (tailnode->A_node->block_loaded_in_GPU == 1)
    {
      tailnode->A_node->block_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
      cublasFree(tailnode->A_node->device_pointer);
#endif
    }

    if (tailnode->B_node->block_loaded_in_GPU == 1)
    {
      tailnode->B_node->block_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
      cublasFree(tailnode->B_node->device_pointer);
#endif
    }

    if (tailnode->C_node->block_loaded_in_GPU == 1)
    {
      tailnode->C_node->block_loaded_in_GPU = 0;
#ifdef HAVE_CUDA
      cublasGetMatrix(tailnode->C_node->M_block, tailnode->C_node->N_block, sizeof(float_t), (void*) tailnode->C_node->device_pointer,
          tailnode->C_node->M_block, (void*) tailnode->C_node->block_dense, tailnode->C_node->M_block);
      cublasFree(tailnode->C_node->device_pointer);
#endif
    }

    tailnode = tailnode->next;
  }

  /* Print statistics. */
  spamm_log("blocks loaded: stream length = %3u cache = %3u A = %3u B = %3u C = %3u total = %4u\n", __FILE__, __LINE__,
      multiply_stream->number_elements, cache_length,
      number_A_blocks_loaded, number_B_blocks_loaded, number_C_blocks_loaded,
      number_A_blocks_loaded+number_B_blocks_loaded+number_C_blocks_loaded);
}

/** Computes the product
 *
 * C = alpha*A*B + beta*C
 *
 * @param alpha the scalar factor multiplying A*B.
 * @param A the matrix A.
 * @param B the matrix B.
 * @param beta the scalar factor multiplying C.
 * @param C the matrix C.
 */
void
spamm_multiply (const float_t alpha, const struct spamm_t *A, const struct spamm_t *B, const float_t beta, struct spamm_t *C)
{
  struct timeval start, stop;
  struct spamm_multiply_stream_t multiply_stream;

  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

  if (A->N != B->M)
  {
    spamm_log("matrix size mismatch, A->N = %i, B->M = %i\n", __FILE__, __LINE__, A->N, B->M);
    exit(1);
  }

  if (A->M != C->M)
  {
    spamm_log("matrix size mismatch, A->M = %i, C->M = %i\n", __FILE__, __LINE__, A->M, C->M);
    exit(1);
  }

  if (B->N != C->N)
  {
    spamm_log("matrix size mismatch, B->N = %i, C->N = %i\n", __FILE__, __LINE__, B->N, C->N);
    exit(1);
  }

  if (A->N_child != B->M_child)
  {
    spamm_log("matrix child size mismatch, A->N_child = %i, B_child->M = %i\n", __FILE__, __LINE__, A->N_child, B->M_child);
    exit(1);
  }

  if (A->M_child != C->M_child)
  {
    spamm_log("matrix child size mismatch, A->M_child = %i, C->M_child = %i\n", __FILE__, __LINE__, A->M_child, C->M_child);
    exit(1);
  }

  if (B->N_child != C->N_child)
  {
    spamm_log("matrix child size mismatch, B->N_child = %i, C->N_child = %i\n", __FILE__, __LINE__, B->N_child, C->N_child);
    exit(1);
  }

  if (A->N_block != B->M_block)
  {
    spamm_log("matrix block size mismatch, A->N_block = %i, B_block->M = %i\n", __FILE__, __LINE__, A->N_block, B->M_block);
    exit(1);
  }

  if (A->M_block != C->M_block)
  {
    spamm_log("matrix block size mismatch, A->M_block = %i, C->M_block = %i\n", __FILE__, __LINE__, A->M_block, C->M_block);
    exit(1);
  }

  if (B->N_block != C->N_block)
  {
    spamm_log("matrix block size mismatch, B->N_block = %i, C->N_block = %i\n", __FILE__, __LINE__, B->N_block, C->N_block);
    exit(1);
  }

  if (C->root != NULL)
  {
    spamm_log("[FIXME] can not handle pre-existing C\n", __FILE__, __LINE__);
    exit(1);
  }

  if (beta != 1.0)
  {
    spamm_log("[FIXME] can not handle (beta = %e) != 1.0\n", __FILE__, __LINE__, beta);
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

  spamm_ll_new(&multiply_stream);

  gettimeofday(&start, NULL);
  spamm_multiply_node(alpha, A->root, B->root, beta, &(C->root), &multiply_stream);
  gettimeofday(&stop, NULL);
  spamm_log("tree recursion: %f s\n", __FILE__, __LINE__, (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

  gettimeofday(&start, NULL);
  spamm_multiply_stream(30000, &multiply_stream);
  gettimeofday(&stop, NULL);
  spamm_log("stream multiply: %f s\n", __FILE__, __LINE__, (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

  spamm_ll_delete(&multiply_stream);
}
