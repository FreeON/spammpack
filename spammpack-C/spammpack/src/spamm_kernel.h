/** @file */

#ifndef __SPAMM_KERNEL_H
#define __SPAMM_KERNEL_H

/** Available stream kernels. */
#define SPAMM_NUMBER_KERNELS 18

/** The different stream kernels.
 */
enum spamm_kernel_t
{
  /** The external sgemm version of the stream kernel. */
  kernel_external_sgemm,

  /** The external sgemm version of the stream kernel with NULL as sgemm_(). */
  kernel_stream_NULL,

  /** The standard stream kernel (SSE). */
  kernel_standard_SSE,

  /** The standard stream kernel, sparse version (SSE). */
  kernel_standard_SSE_sparse,

  /** The standard stream kernel without norm checks (SSE). */
  kernel_standard_no_checks_SSE,

  /** The standard stream kernel (SSE4.1). */
  kernel_standard_SSE4_1,

  /** The standard stream kernel without norm checks (SSE4.1). */
  kernel_standard_no_checks_SSE4_1,

  /** The Z-curve stream kernel (SSE). */
  kernel_Z_curve_SSE,

  /** The Z-curve stream kernel without norm checks (SSE). */
  kernel_Z_curve_no_checks_SSE,

  /** The Z-curve stream kernel (SSE4.1). */
  kernel_Z_curve_SSE4_1,

  /** The Z-curve stream kernel without norm checks (SSE4.1). */
  kernel_Z_curve_no_checks_SSE4_1,

  /** The hierarchical stream kernel (SSE). */
  kernel_hierarchical_SSE,

  /** The hierarchical stream kernel without norm checks (SSE). */
  kernel_hierarchical_no_checks_SSE,

  /** The hierarchical stream kernel (SSE4.1). */
  kernel_hierarchical_SSE4_1,

  /** The hierarchical stream kernel without norm checks (SSE4.1). */
  kernel_hierarchical_no_checks_SSE4_1,

  /** The hierarchical Z-curve stream kernel (SSE). */
  kernel_hierarchical_Z_curve_SSE,

  /** The hierarchical Z-curve stream kernel without norm checks (SSE). */
  kernel_hierarchical_Z_curve_no_checks_SSE,

  /** The hierarchical Z-curve stream kernel (SSE4.1). */
  kernel_hierarchical_Z_curve_SSE4_1,

  /** The hierarchical Z-curve stream kernel without norm checks (SSE4.1). */
  kernel_hierarchical_Z_curve_no_checks_SSE4_1
};

/** The basic information in a stream element.
 */
struct spamm_multiply_stream_t
{
  /** A pointer to the kernel tier matrix node of A. */
  struct spamm_data_t *A;

  /** A pointer to the kernel tier matrix node of B. */
  struct spamm_data_t *B;

  /** A pointer to the kernel tier matrix node of C. */
  struct spamm_data_t *C;
};

const char *
spamm_kernel_get_name (const unsigned int i);

enum spamm_kernel_t
spamm_kernel_get_kernel (const char* name);

enum spamm_layout_t
spamm_kernel_suggest_layout (const enum spamm_kernel_t kernel);

enum spamm_layout_t
spamm_kernel_get_layout (const char *name);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_SSE (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_SSE_sparse (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream without norm checks.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_no_checks_SSE (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_SSE4_1 (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream without norm checks.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_no_checks_SSE4_1 (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_Z_curve_SSE (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_Z_curve_no_checks_SSE (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_Z_curve_SSE4_1 (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_Z_curve_no_checks_SSE4_1 (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_hierarchical_SSE (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_hierarchical_no_checks_SSE (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_hierarchical_SSE4_1 (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_hierarchical_no_checks_SSE4_1 (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_hierarchical_Z_curve_SSE (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_hierarchical_Z_curve_no_checks_SSE (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_hierarchical_Z_curve_SSE4_1 (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel_hierarchical_Z_curve_no_checks_SSE4_1 (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

/** Process the multiply stream. Version written in C.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_external_sgemm (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream,
    const short call_external_sgemm);

#endif
