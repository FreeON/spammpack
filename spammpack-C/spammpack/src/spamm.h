/** @file */

#if ! defined(__SPAMM_H)

/** Define in case spamm.h has been included. */
#define __SPAMM_H 1

#include "config.h"
#include "spamm_config.h"
#include "spamm_kernel.h"
#include <stdio.h>

/** Definition of global floating point precision.
 *
 * This is either float or double, as defined by the configure script.
 */
typedef FLOATING_PRECISION floating_point_t;

/** The severity levels for the logger.
 */
enum spamm_log_severity_t
{
  /** A fatal message that should always be displayed. */
  fatal,

  /** A message that might be printed if the global loglevel is set to at
   * least info. */
  info,

  /** A message that only gets printed if the global loglevel is set to debug.
   */
  debug
};

/** Define shortcut macro for logging.
 *
 * Typical use of this macro:
 *
 * <code>LOG_FATAL("opening new file: %s\n", filename);</code>
 *
 * Compare to spamm_log().
 *
 * @param format Format string. See printf() for a detailed description of its
 *               syntax.
 */
#define LOG_FATAL(format, ...) spamm_log(fatal, format, __FILE__, __LINE__, __VA_ARGS__)

/** Define shortcut macro for logging.
 *
 * Typical use of this macro:
 *
 * <code>LOG2_FATAL("a message without arguments\n");</code>
 *
 * Compare to spamm_log().
 *
 * @param format Format string. See printf() for a detailed description of its
 *               syntax.
 */
#define LOG2_FATAL(format) spamm_log(fatal, format, __FILE__, __LINE__)

/** Define shortcut macro for logging.
 *
 * Typical use of this macro:
 *
 * <code>LOG_INFO("opening new file: %s\n", filename);</code>
 *
 * Compare to spamm_log().
 *
 * @param format Format string. See printf() for a detailed description of its
 *               syntax.
 */
#ifdef ENABLE_LOGGING
#define LOG_INFO(format, ...) spamm_log(info, format, __FILE__, __LINE__, __VA_ARGS__)
#else
#define LOG_INFO(format, ...)
#endif

/** Define shortcut macro for logging.
 *
 * Typical use of this macro:
 *
 * <code>LOG2_INFO("a message without arguments\n");</code>
 *
 * Compare to spamm_log().
 *
 * @param format Format string. See printf() for a detailed description of its
 *               syntax.
 */
#ifdef ENABLE_LOGGING
#define LOG2_INFO(format) spamm_log(info, format, __FILE__, __LINE__)
#else
#define LOG2_INFO(format)
#endif

/** Define shortcut macro for logging.
 *
 * Typical use of this macro:
 *
 * <code>LOG_DEBUG("opening new file: %s\n", filename);</code>
 *
 * Compare to spamm_log().
 *
 * @param format Format string. See printf() for a detailed description of its
 *               syntax.
 */
#ifdef ENABLE_LOGGING
#define LOG_DEBUG(format, ...) spamm_log(debug, format, __FILE__, __LINE__, __VA_ARGS__)
#else
#define LOG_DEBUG(format, ...)
#endif

/** Define shortcut macro for logging.
 *
 * Typical use of this macro:
 *
 * <code>LOG2_DEBUG("a message without arguments\n");</code>
 *
 * Compare to spamm_log().
 *
 * @param format Format string. See printf() for a detailed description of its
 *               syntax.
 */
#ifdef ENABLE_LOGGING
#define LOG2_DEBUG(format) spamm_log(debug, format, __FILE__, __LINE__)
#else
#define LOG2_DEBUG(format)
#endif

/* Definition of return codes. */

/** Return code: Everything went fine. */
#define SPAMM_RESULT_OK 0

/** Return code: Something went wrong. */
#define SPAMM_RESULT_FAILED -1

/** Return code: The result was below the threshold. */
#define SPAMM_RESULT_BELOW_THRESHOLD 1

/** Return code: Trying to set a matrix element which is zero. */
#define SPAMM_RESULT_ZERO_ELEMENT 2

/** The basic matrix data type.
 */
struct spamm_t
{
  /** Number of rows of matrix. */
  unsigned int M;

  /** Number of columns of matrix. */
  unsigned int N;

  /** Padded size of matrix.
   *
   * Since we are using a square matrix for the basic matrix block and a
   * square matrix for the children node matrix, the padded matrix can only be
   * square.
   */
  unsigned int N_padded;

  /** The Frobenius norm of the matrix underneath this node. */
  floating_point_t norm;

  /** The depth of the tree. */
  unsigned int tree_depth;

  /** The kernel tier determines at which point the SpAMM kernel will get
   * involved.
   *
   * At and below the kernel tier, the dense matrix blocks are allocated
   * contiguously.
   */
  unsigned int kernel_tier;

  /** The number of non-zero blocks. */
  unsigned int number_nonzero_blocks;

  /** The root node. */
  struct spamm_node_t *root;
};

/** A node in the tree.
 *
 * This structure describes a node in the tree.
 */
struct spamm_node_t
{
  /** The tree tier.
   *
   * The root is tier 0, the last tier is equal to the tree_depth.
   */
  unsigned int tier;

  /** The depth of the tree. */
  unsigned int tree_depth;

  /** The kernel tier determines at which point the SpAMM kernel will get
   * involved.
   *
   * At and below the kernel tier, the tree is allocated contiguously.
   */
  unsigned int kernel_tier;

  /** The rows of the padded matrix covered in this node.
   *
   * The indices are meant as [M_lower, M_upper[, i.e. the upper limit is not
   * included in the interval.
   */
  unsigned int M_lower;

  /** The rows of the padded matrix covered in this node.
   *
   * The indices are meant as [M_lower, M_upper[, i.e. the upper limit is not
   * included in the interval.
   */
  unsigned int M_upper;

  /** The columns of the padded matrix covered in this node.
   *
   * The indices are meant as [N_lower, N_upper[, i.e. the upper limit is not
   * included in the interval.
   */
  unsigned int N_lower;

  /** The columns of the padded matrix covered in this node.
   *
   * The indices are meant as [N_lower, N_upper[, i.e. the upper limit is not
   * included in the interval.
   */
  unsigned int N_upper;

  unsigned int M_lower_kernel_tier;
  unsigned int M_upper_kernel_tier;
  unsigned int N_lower_kernel_tier;
  unsigned int N_upper_kernel_tier;

  /** The linear index of this block along the curve. */
  unsigned int index;

  /** The Frobenius norm of the matrix underneath this node. */
  floating_point_t norm;

  /** The square of the Frobenius norm of the matrix underneath this node. */
  floating_point_t norm2;

  /** At the non-block level, pointers to the children nodes.
   *
   * The pointers can be accessed in 2 ways:
   *
   * - As a 2-D array using spamm_dense_index()
   * - As a 1-D array which is sorted on the index.
   */
  struct spamm_node_t *child[SPAMM_N_CHILD][SPAMM_N_CHILD];

  /** At the block level, the dense matrix data. */
  floating_point_t *block_dense;

  /** At the block level, the dilated dense matrix data. */
  floating_point_t *block_dense_dilated;
};

/** Linear quadtree.
 */
struct spamm_linear_quadtree_t
{
  /** The linear quadtree index of this data block.
   */
  unsigned int index;

  /** The number of rows of the dense matrix block. */
  unsigned int M;

  /** The number of columns of the dense data block. */
  unsigned int N;

  /** The data.
   */
  floating_point_t *block_dense;
};

/** Tree statistics.
 *
 * This structure is the result of a call to spamm_tree_stats().
 */
struct spamm_tree_stats_t
{
  /** The number of nodes. */
  unsigned int number_nodes;

  /** The number of dense blocks. */
  unsigned int number_dense_blocks;

  /** The memory consumption of the tree. */
  unsigned int memory_tree;

  /** The memory consumption of the dense blocks. */
  unsigned int memory_dense_blocks;

  /** The memory consumption total. */
  unsigned int memory_total;

  /** The average sparsity of the dense blocks. */
  floating_point_t average_sparsity;
};

/* Function declarations. */

void
spamm_add (const floating_point_t alpha, const struct spamm_t *A, const floating_point_t beta, struct spamm_t *B);

void
spamm_add_node (const floating_point_t alpha, const struct spamm_node_t *A_node, const floating_point_t beta, struct spamm_node_t **B_node);

void *
spamm_allocate (size_t size);

unsigned int
spamm_block_index (const unsigned int i, const unsigned int j,
    const unsigned int M_block, const unsigned int N_block,
    const unsigned int M_kernel, const unsigned int N_kernel);

void
spamm_delete (struct spamm_t *A);

void
spamm_delete_node (struct spamm_node_t **node);

int
spamm_dense_index (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N);

void
spamm_dense_to_spamm (const unsigned int M, const unsigned int N,
    const char operation,
    const floating_point_t *A_dense, struct spamm_t *A);

void
spamm_free (void *data);

floating_point_t
spamm_get (const unsigned int i, const unsigned int j, const struct spamm_t *A);

enum spamm_log_severity_t
spamm_get_loglevel ();

void
spamm_int_to_binary (const unsigned int integer, const int width, char *binary_string);

unsigned int
spamm_linear_index_2D (const unsigned int i, const unsigned int j);

void
spamm_linear_to_coordinates (const unsigned int index, unsigned int *i,
    unsigned int *j, const unsigned int M, const unsigned int N,
    const unsigned int M_block, const unsigned int N_block);

void
spamm_log (const enum spamm_log_severity_t severity, const char *format,
    const char *filename, const unsigned int linenumber, ...);

unsigned int
spamm_multiply (const floating_point_t tolerance,
    const floating_point_t alpha, const struct spamm_t *A,
    const struct spamm_t *B, const floating_point_t beta, struct spamm_t *C);

void
spamm_multiply_scalar (const floating_point_t alpha, struct spamm_t *A);

void
spamm_new (const unsigned int M, const unsigned int N, struct spamm_t *A);

struct spamm_node_t *
spamm_new_childnode (const unsigned int tier,
    const unsigned int tree_depth,
    const unsigned int M_lower, const unsigned int M_upper,
    const unsigned int N_lower, const unsigned int N_upper,
    const unsigned int M_lower_kernel_tier, const unsigned int M_upper_kernel_tier,
    const unsigned int N_lower_kernel_tier, const unsigned int N_upper_kernel_tier,
    const unsigned int kernel_tier,
    floating_point_t *block_dense, floating_point_t *block_dense_dilated);

struct spamm_node_t *
spamm_new_node ();

unsigned int
spamm_number_nonzero (const struct spamm_t *A);

void
spamm_print_dense (const unsigned int M, const unsigned int N, const floating_point_t *A_dense);

void
spamm_print_node (const struct spamm_node_t *node);

void
spamm_print_spamm (const struct spamm_t *A);

void
spamm_print_tree (const struct spamm_t *A);

void
spamm_read_MM (const char *filename, struct spamm_t *A);

int
spamm_set (const unsigned int i, const unsigned int j, const floating_point_t Aij, struct spamm_t *A);

void
spamm_set_loglevel (const enum spamm_log_severity_t loglevel);

void
spamm_sgemm_trivial (const char opA, const char opB,
    const unsigned int M, const unsigned int N, const unsigned int K,
    const floating_point_t alpha,
    const floating_point_t *A_block_dense,
    const unsigned int lda,
    const floating_point_t *B_block_dense,
    const unsigned int ldb,
    const floating_point_t beta,
    floating_point_t *C_block_dense,
    const unsigned int ldc);

void
spamm_spamm_to_dense (const struct spamm_t *A, floating_point_t **A_dense);

void
spamm_tree_stats (struct spamm_tree_stats_t *stats, const struct spamm_t *A);

char *
spamm_version ();

#endif
