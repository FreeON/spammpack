/** @file */

#ifndef __SPAMM_KERNEL_H
#define __SPAMM_KERNEL_H

/** The basic information in a stream element.
 */
struct multiply_stream_t
{
  float *A_block;
  float *B_block;
  float *C_block;
  float  norm[32];
};

void
spamm_stream_kernel (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct multiply_stream_t *multiply_stream);

#endif
