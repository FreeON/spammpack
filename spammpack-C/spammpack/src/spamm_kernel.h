/** @file */

#ifndef __SPAMM_KERNEL_H
#define __SPAMM_KERNEL_H

/** The basic information in a stream element.
 */
struct spamm_multiply_stream_t
{
  /** A pointer to the kernel dense matrix block of A. This block is assumed
   * to be dilated by 4 for SSE processing. */
  float *A_block;

  /** A pointer to the kernel dense matrix block of B. */
  float *B_block;

  /** A pointer to the kernel dense matrix block of B. */
  float *C_block;

  /** The 2x16 norms of the underlying matrix blocks of A and B. */
  float  norm[32];
};

/** Process the multiply stream.
 *
 * @param number_stream_elements The size of the multiply stream.
 * @param alpha The factor \f$\alpha\f$.
 * @param tolerance The SpAMM tolerance.
 * @param multiply_stream The multiply stream.
 */
void
spamm_stream_kernel (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct spamm_multiply_stream_t *multiply_stream);

#endif
