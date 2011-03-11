/** @file */

#ifndef __SPAMM_KERNEL_H
#define __SPAMM_KERNEL_H

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
spamm_stream_kernel_no_checks_SSE (const unsigned int number_stream_elements,
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
spamm_stream_kernel_C (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

#endif
