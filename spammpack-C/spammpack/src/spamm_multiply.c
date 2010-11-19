#include "spamm.h"
#include "config.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

enum hll_order
{
  row_major, column_major
};

/** \private Creates the horizontal link layer.
 *
 * @param order The order of the hll.
 * @param A The matrix tree.
 *
 * @return The number of blocks in the hll.
 */
void
spamm_multiply_create_hll (const enum hll_order order, struct spamm_t *A, unsigned int *number_nodes, struct spamm_node_t **first_node)
{
  unsigned int i, j;
  struct spamm_node_t *node;
  struct spamm_node_t *last_node = NULL;

  /* Reset first node. */
  *first_node = NULL;
  *number_nodes = 0;

  /* Construct hll. */
  switch (order)
  {
    case row_major:
      for (i = 0; i < A->M; i += SPAMM_N_KERNEL) {
        for (j = 0; j < A->N; j += SPAMM_N_KERNEL)
        {
          //LOG_DEBUG("finding node of A(%u,%u)\n", i, j);
          node = spamm_get_node(i, j, A->kernel_tier, A);

          if (node != NULL)
          {
            //spamm_print_node(node);
            (*number_nodes)++;

            /* Link to previous node. */
            if (last_node != NULL)
            {
              last_node->next = node;
            }

            else
            {
              *first_node = node;
            }

            last_node = node;
          }
        }
      }
      break;

    case column_major:
      for (j = 0; j < A->N; j += SPAMM_N_KERNEL) {
        for (i = 0; i < A->M; i += SPAMM_N_KERNEL)
        {
          //LOG_DEBUG("finding node of A(%u,%u)\n", i, j);
          node = spamm_get_node(i, j, A->kernel_tier, A);

          if (node != NULL)
          {
            //spamm_print_node(node);
            (*number_nodes)++;

            /* Link to previous node. */
            if (last_node != NULL)
            {
              last_node->next = node;
            }

            else
            {
              *first_node = node;
            }

            last_node = node;
          }
        }
      }
      break;
  }
}

/** \private Generate convolution.
 */
unsigned int
spamm_multiply_convolute (const unsigned int number_hll_A,
    const unsigned int *hll_nodelist_A_index,
    floating_point_t **hll_nodelist_A_block,
    const unsigned int number_hll_B,
    const unsigned int *hll_nodelist_B_index,
    floating_point_t **hll_nodelist_B_block,
    const unsigned int max_number_elements,
    struct multiply_stream_t *multiply_stream,
    unsigned int *multiply_stream_index)
{
  unsigned int i, j, k;
  unsigned int stream_index;

  /* [FIXME] This step should be made O(N^2). We need to introduce some
   * linking in the the nodelists, sort of like a skip-list, so we can loop
   * over ranges of fixed k.
   */
  stream_index = 0;
  for (i = 0; i < number_hll_A; i++) {
    for (j = 0; j < number_hll_B; j++)
    {
      if ((hll_nodelist_A_index[i] & 0x12492492) == (hll_nodelist_B_index[j] & 0x12492492))
      {
        multiply_stream_index[stream_index] = hll_nodelist_A_index[i] | hll_nodelist_B_index[j];
        multiply_stream[stream_index].A_block = hll_nodelist_A_block[i];
        multiply_stream[stream_index].B_block = hll_nodelist_B_block[j];

        /* [FIXME] The norms should be loaded from the actual matrix blocks.
         * Format: A(1,1), A(1,2), A(1,3), A(1,4), A(2,1), ...., B(1,1),
         * B(1,2) ....
         */
        for (k = 0; k < 32; k++)
        {
          multiply_stream[stream_index].norm[k] = 1.0;
        }

        stream_index++;
        if (stream_index > max_number_elements)
        {
          LOG2_FATAL("insufficient space in multiply_stream\n");
          exit(1);
        }
      }
    }
  }

  return stream_index;
}

/** \private Build stream.
 *
 * The stream is built by allocating the proper nodes in the C tree and
 * linking the dense blocks in C into the stream.
 */
void
spamm_multiply_build_stream (const unsigned int number_elements,
    struct multiply_stream_t *multiply_stream,
    unsigned int *multiply_stream_index,
    struct spamm_t *C)
{
  short bit_index, masked_Cij;
  unsigned int index;
  unsigned int Cij;
  unsigned int i, j;
  char index_string[1000];
  char Cij_string[1000];

  /* Construct indices of C blocks, allocate C tree, and assign C dense blocks
   * to multiply stream.
   */
  for (index = 0; index < number_elements; index++)
  {
    Cij = multiply_stream_index[index] & 0x2db6db6d;

    for (i = 0, j = 0, bit_index = 0; bit_index < 10; bit_index++)
    {
      masked_Cij = Cij & 7;
      switch (masked_Cij)
      {
        case 1:
          j |= 1 << bit_index;
          break;

        case 4:
          i |= 1 << bit_index;
          break;

        case 5:
          i |= 1 << bit_index;
          j |= 1 << bit_index;
          break;
      }
      Cij >>= 3;
    }

    i *= SPAMM_N_BLOCK;
    j *= SPAMM_N_BLOCK;

    //spamm_int_to_binary(multiply_stream_index[index], 32, index_string);
    //spamm_int_to_binary(multiply_stream_index[index] & 0x2db6db6d, 32, Cij_string);
    //LOG_INFO("analyzing stream element %06u: %s --> %s --> (%u,%u)\n", index, index_string, Cij_string, i, j);

    /* Create necessary nodes in C and link into stream. */
    spamm_create_tree(i, j, C);
  }
}


/** \private Computes
 *
 * C_node = alpha*A_node*B_node + C_node
 *
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
spamm_multiply_node (const floating_point_t tolerance,
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

  if (A_node == NULL || B_node == NULL)
  {
    /* Nothing to do here. */
    return number_products;
  }

  if (A_node->tier != B_node->tier)
  {
    LOG_FATAL("There seems to be a bug here, A_node-tier != B_node->tier (%u != %u)\n", A_node->tier, B_node->tier);
    exit(1);
  }

  /* Create new node. */
  if (*C_node == NULL)
  {
    if (A_node->tier == 0)
    {
      LOG2_DEBUG("creating new C root node\n");

      *C_node = spamm_new_childnode(0, A_node->tree_depth,
          A_node->M_lower, A_node->N_upper, B_node->N_lower, B_node->N_upper,
          0, 0, 0, 0,
          A_node->kernel_tier, NULL, NULL);
    }

    else
    {
      LOG2_FATAL("C_node is not set but we are not at the root\n");
      exit(1);
    }
  }

  /* Do some work if at the kernel tier. */
  if (A_node->tier == A_node->kernel_tier)
  {
    /* Pop this product onto multiply stream. */
    multiply_stream[*number_multiply_stream_elements].A_block = A_node->block_dense_dilated;
    multiply_stream[*number_multiply_stream_elements].B_block = B_node->block_dense;
    multiply_stream[*number_multiply_stream_elements].C_block = (*C_node)->block_dense;
    for (l = 0; l < 32; l++)
    {
      multiply_stream[*number_multiply_stream_elements].norm[l] = 1;
    }
    (*number_multiply_stream_elements)++;
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
              (*C_node)->child[i][j] = spamm_new_childnode((*C_node)->tier+1, (*C_node)->tree_depth,
                  (*C_node)->M_lower+i*((*C_node)->M_upper-(*C_node)->M_lower)/SPAMM_N_CHILD,
                  (*C_node)->M_lower+(i+1)*((*C_node)->M_upper-(*C_node)->M_lower)/SPAMM_N_CHILD,
                  (*C_node)->N_lower+j*((*C_node)->N_upper-(*C_node)->N_lower)/SPAMM_N_CHILD,
                  (*C_node)->N_lower+(j+1)*((*C_node)->N_upper-(*C_node)->N_lower)/SPAMM_N_CHILD,
                  (*C_node)->M_lower_kernel_tier, (*C_node)->M_upper_kernel_tier,
                  (*C_node)->N_lower_kernel_tier, (*C_node)->N_upper_kernel_tier,
                  (*C_node)->kernel_tier, (*C_node)->block_dense, (*C_node)->block_dense_dilated);
            }

            /* Recurse. */
            number_products += spamm_multiply_node(tolerance, alpha,
                A_node->child[i][k], B_node->child[k][j], &(*C_node)->child[i][j],
                number_multiply_stream_elements, multiply_stream);
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

/** Computes the product
 *
 * \f$C = \alpha A \times B + \beta C\f$
 *
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
spamm_multiply (floating_point_t tolerance,
    const floating_point_t alpha, struct spamm_t *A,
    struct spamm_t *B, const floating_point_t beta, struct spamm_t *C)
{
  struct timeval hll_start, hll_stop;
  struct timeval tree_start, tree_stop;
  struct timeval convolution_start, convolution_stop;
  struct timeval streambuild_start, streambuild_stop;
  struct timeval stream_start, stream_stop;
  struct timeval total_start, total_stop;

  struct multiply_stream_t *multiply_stream;
  unsigned int *multiply_stream_index;

  unsigned int *hll_nodelist_A_index;
  unsigned int *hll_nodelist_B_index;
  floating_point_t **hll_nodelist_A_block;
  floating_point_t **hll_nodelist_B_block;

  struct spamm_node_t *first_hll_node_A;
  struct spamm_node_t *first_hll_node_B;
  struct spamm_node_t *node;

  double hll_time;
  double tree_time;
  double convolution_time;
  double streambuild_time;
  double stream_time;

  unsigned int number_multiply_stream_elements = 0;
  unsigned int number_products = 0;
  unsigned int number_hll_A, number_hll_B;
  unsigned int hll_index;

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
  multiply_stream_index = (unsigned int*) malloc(sizeof(unsigned int)*max_number_stream_elements);

  /* Set up horizontal link list layer. */
  gettimeofday(&hll_start, NULL);
  spamm_multiply_create_hll(column_major, A, &number_hll_A, &first_hll_node_A);
  spamm_multiply_create_hll(row_major, B, &number_hll_B, &first_hll_node_B);
  hll_nodelist_A_index = (unsigned int*) malloc(sizeof(unsigned int)*number_hll_A);
  hll_nodelist_A_block = (floating_point_t**) malloc(sizeof(floating_point_t*)*number_hll_A);
  for (node = first_hll_node_A, hll_index = 0; node != NULL; node = node->next)
  {
    hll_nodelist_A_index[hll_index] = node->index_3D_column;
    hll_nodelist_A_block[hll_index++] = node->block_dense_dilated;
  }
  hll_nodelist_B_index = (unsigned int*) malloc(sizeof(unsigned int)*number_hll_B);
  hll_nodelist_B_block = (floating_point_t**) malloc(sizeof(floating_point_t*)*number_hll_B);
  for (node = first_hll_node_B, hll_index = 0; node != NULL; node = node->next)
  {
    hll_nodelist_B_index[hll_index] = node->index_3D_row;
    hll_nodelist_B_block[hll_index++] = node->block_dense;
  }
  gettimeofday(&hll_stop, NULL);
  hll_time = (hll_stop.tv_sec-hll_start.tv_sec)+(hll_stop.tv_usec-hll_start.tv_usec)/(double) 1e6;
  printf("hll setup: hll A = %u elements, hll B = %u elements, %f s\n", number_hll_A, number_hll_B, hll_time);

  /* Construct convolution. */
  gettimeofday(&convolution_start, NULL);
  number_products = spamm_multiply_convolute(number_hll_A, hll_nodelist_A_index, hll_nodelist_A_block,
      number_hll_B, hll_nodelist_B_index, hll_nodelist_B_block, max_number_stream_elements, multiply_stream, multiply_stream_index);
  gettimeofday(&convolution_stop, NULL);
  convolution_time = (convolution_stop.tv_sec-convolution_start.tv_sec)+(convolution_stop.tv_usec-convolution_start.tv_usec)/(double) 1e6;
  printf("convolution: constructed %u elements, %f s\n", number_products, convolution_time);

  /* Link in C tree. */
  gettimeofday(&streambuild_start, NULL);
  spamm_multiply_build_stream(number_products, multiply_stream, multiply_stream_index, C);
  gettimeofday(&streambuild_stop, NULL);
  streambuild_time = (streambuild_stop.tv_sec-streambuild_start.tv_sec)+(streambuild_stop.tv_usec-streambuild_start.tv_usec)/(double) 1e6;
  printf("streambuild: %f s\n", streambuild_time);

  //gettimeofday(&tree_start, NULL);
  //number_products = spamm_multiply_node(tolerance, alpha, A->root, B->root, &(C->root), &number_multiply_stream_elements, multiply_stream);
  //gettimeofday(&tree_stop, NULL);
  //tree_time = (tree_stop.tv_sec-tree_start.tv_sec)+(tree_stop.tv_usec-tree_start.tv_usec)/(double) 1e6;
  //printf("symbolic multiply: placed %u elements into stream, %f s\n", number_multiply_stream_elements, tree_time);

  gettimeofday(&stream_start, NULL);
  spamm_stream_kernel(number_multiply_stream_elements, alpha, tolerance, multiply_stream);
  gettimeofday(&stream_stop, NULL);
  stream_time = (stream_stop.tv_sec-stream_start.tv_sec)+(stream_stop.tv_usec-stream_start.tv_usec)/(double) 1e6;
  printf("stream multiply: %f s\n", stream_time);

  printf("symbolic is %1.1f times stream\n", tree_time/stream_time);

  free(multiply_stream);

  gettimeofday(&total_stop, NULL);
  printf("total time for multiply: %f s\n", (total_stop.tv_sec-total_start.tv_sec)+(total_stop.tv_usec-total_start.tv_usec)/(double) 1e6);

  return number_products;
}
