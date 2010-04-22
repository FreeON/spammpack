#include "config.h"
#include <stdio.h>

typedef FLOATING_PRECISION float_t;

struct spamm_t
{
  /* Matrix dimensions. */
  int M, N;

  /* The padded dimensions. */
  int M_padded, N_padded;

  /* The dimensions of the data blocks at the leaf level. */
  int M_block, N_block;

  /* The dimensions of subdivisions on each child node. */
  int M_child, N_child;

  /* The matrix element threshold. Elements below this threshold are not
   * stored. */
  float_t threshold;

  /* The root node. */
  struct spamm_node_t *root;
};

enum spamm_block_ordering_t
{
  P, Q, R, S
};

struct spamm_node_t
{
  /* The dimensions covered in the padded matrix in this node. The indices are
   * meant as [M_lower, M_upper[, i.e. the upper limit is not included in the
   * interval. */
  int M_lower, M_upper;
  int N_lower, N_upper;

  /* The dimensions of the data blocks at the leaf level. */
  int M_block, N_block;

  /* The dimensions of subdivisions on each child node. */
  int M_child, N_child;

  /* The matrix element threshold. Elements below this threshold are not
   * stored. */
  float_t threshold;

  /* The linear index of this block along the curve. */
  unsigned int index;

  /* The name of the block ordering pattern. */
  enum spamm_block_ordering_t ordering;

  /* At the non-block level, pointers to the children nodes. */
  struct spamm_node_t **child;

  /* At the block level, the dense matrix data. */
  float_t *block_dense;
};

struct spamm_tree_stats_t
{
  /* The number of nodes. */
  int number_nodes;

  /* The number of dense blocks. */
  int number_dense_blocks;

  /* The memory consumption of the tree. */
  int memory_tree;

  /* The memory consumption of the dense blocks. */
  int memory_dense_blocks;
};

struct spamm_multiply_stream_t
{
  /* Links to the first and last node in list. */
  struct spamm_multiply_stream_node_t *first;
  struct spamm_multiply_stream_node_t *last;
};

struct spamm_multiply_stream_node_t
{
  /* Links to the previous and the next element. */
  struct spamm_multiply_stream_node_t *previous;
  struct spamm_multiply_stream_node_t *next;

  /* Index triple of matrices. */
  unsigned int A_index;
  unsigned int B_index;
  unsigned int C_index;

  /* Pointers into nodes corresponding to indices. */
  const struct spamm_node_t *A_node;
  const struct spamm_node_t *B_node;
  const struct spamm_node_t *C_node;
};

void
spamm_log (const char *format, const char *filename, const int linenumber, ...);

int
spamm_dense_index (const int i, const int j, const int M, const int N);

void
spamm_new (const int M, const int N, const int M_block, const int N_block,
    const int M_child, const int N_child, const float_t threshold,
    struct spamm_t *A);

void
spamm_new_node (struct spamm_node_t **node);

void
spamm_delete (struct spamm_t *A);

void
spamm_dense_to_spamm (const int M, const int N, const int M_block,
    const int N_block, const int M_child, const int N_child,
    const float_t threshold, const float_t *A_dense, struct spamm_t *A);

void
spamm_spamm_to_dense (const struct spamm_t *A, float_t **A_dense);

float_t
spamm_get (const int i, const int j, const struct spamm_t *A);

void
spamm_set (const int i, const int j, const float_t Aij, struct spamm_t *A);

void
spamm_print_dense (const int M, const int N, const float_t *A_dense);

void
spamm_print_spamm (const struct spamm_t *A);

void
spamm_print_node (const struct spamm_node_t *node);

void
spamm_print_tree (const struct spamm_t *A);

void
spamm_add (const float_t alpha, const struct spamm_t *A, const float_t beta, struct spamm_t *B);

void
spamm_multiply (const float_t alpha, const struct spamm_t *A, const struct spamm_t *B, const float_t beta, struct spamm_t *C);

void
spamm_tree_stats (struct spamm_tree_stats_t *stats, const struct spamm_t *A);

void
spamm_read_MM (const char *filename, const int M_block, const int N_block,
    const int M_child, const int N_child, const float_t threshold,
    struct spamm_t *A);

void
spamm_ll_new (struct spamm_multiply_stream_t *list);

void
spamm_ll_delete (struct spamm_multiply_stream_t *list);

void
spamm_ll_append (const unsigned int A_index, const struct spamm_node_t *A_node,
    const unsigned int B_index, const struct spamm_node_t *B_node,
    const unsigned int C_index, const struct spamm_node_t *C_node,
    struct spamm_multiply_stream_t *list);

void
spamm_ll_sort (struct spamm_multiply_stream_t *list);

void
spamm_ll_print (struct spamm_multiply_stream_t *list);
