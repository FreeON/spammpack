#include "spamm.h"
#include "config.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

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
    unsigned int *number_multiply_stream_elements,
    struct multiply_stream_t *multiply_stream)
{
  floating_point_t beta = 1.0;
  int i, j, k, l, m;
  char bitstring_A[100];
  char bitstring_B[100];
  char bitstring_C[100];
  unsigned int mask_A, mask_B, mask_C;
  struct spamm_node_t *C_child_node;
  struct spamm_multiply_stream_element_t *multiply_stream_element;
  struct spamm_ll_iterator_t *iterator_A, *iterator_B, *iterator_C;
  struct spamm_ll_node_t *linear_node_A, *linear_node_B, *linear_node_C, *linear_node_C_next;
  struct spamm_linear_quadtree_t *linear_A, *linear_B, *linear_C, *linear_C_next;
  unsigned int number_products = 0;
  unsigned int kernel_block_N;

  assert(A_node != NULL);
  assert(B_node != NULL);
  assert(A_node->tier == B_node->tier);

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

  /* Recurse down the tree. */
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
                C_child_node->block_dense = (floating_point_t*) spamm_allocate(sizeof(floating_point_t)*SPAMM_N_KERNEL*SPAMM_N_KERNEL);
                for (l = 0; l < SPAMM_N_KERNEL; ++l) {
                  for (m = 0; m < SPAMM_N_KERNEL; ++m)
                  {
                    C_child_node->block_dense[spamm_dense_index(l, m, SPAMM_N_KERNEL, SPAMM_N_KERNEL)] = 0.0;
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

            if (A_node->tier < A_node->kernel_tier)
            {
              /* Recurse. */
              number_products += spamm_multiply_node(algorithm, tolerance, alpha,
                  A_node->child[i][k], B_node->child[k][j], &(*C_node)->child[i][j],
                  number_multiply_stream_elements, multiply_stream);
            }

            else if (A_node->tier == A_node->kernel_tier)
            {
              switch (algorithm)
              {
                case tree:
                  break;

                case cache:
                  /* Pop this product onto multiply stream. */
                  multiply_stream[*number_multiply_stream_elements].A_block = A_node->block_dense;
                  multiply_stream[*number_multiply_stream_elements].B_block = B_node->block_dense;
                  multiply_stream[*number_multiply_stream_elements].C_block = (*C_node)->block_dense;
                  for (l = 0; l < 32; l++)
                  {
                    multiply_stream[*number_multiply_stream_elements].norm[l] = 1;
                  }
                  (*number_multiply_stream_elements)++;
                  break;

                default:
                  LOG2_FATAL("[FIXME]\n");
                  exit(1);
                  break;
              }
            }
          }

          else
          {
            LOG2_DEBUG("dropping product, it is below tolerance\n");
          }
        }
      }
    }
  }

  return number_products;
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
  struct multiply_stream_t *multiply_stream;
  unsigned int number_multiply_stream_elements = 0;
  unsigned int number_products = 0;

  unsigned int max_number_stream_elements;
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
      number_products = spamm_multiply_node(algorithm, tolerance, alpha, A->root, B->root, &(C->root), NULL, NULL);
      break;

    case cache:
      LOG2_INFO("using cache algorithm\n");
    case cache_redundant:
      if (algorithm == cache_redundant)
      {
        LOG2_INFO("using cache (redundant) algorithm\n");
      }

      max_number_stream_elements = pow(A->N_padded/SPAMM_N_BLOCK/pow(2, SPAMM_KERNEL_DEPTH), 3);
      max_memory = max_number_stream_elements*sizeof(struct multiply_stream_t);

      LOG_INFO("max number of stream elements: %u\n", max_number_stream_elements);

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
      multiply_stream = (struct multiply_stream_t*) malloc(sizeof(struct multiply_stream_t)*max_number_stream_elements);

      gettimeofday(&start, NULL);
      number_products = spamm_multiply_node(algorithm, tolerance, alpha, A->root, B->root, &(C->root), &number_multiply_stream_elements, multiply_stream);
      gettimeofday(&stop, NULL);
      LOG_INFO("symbolic multiply: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

      gettimeofday(&start, NULL);
      spamm_stream_kernel(number_multiply_stream_elements, alpha, tolerance, multiply_stream);
      gettimeofday(&stop, NULL);
      LOG_INFO("stream multiply: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

      free(multiply_stream);

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
