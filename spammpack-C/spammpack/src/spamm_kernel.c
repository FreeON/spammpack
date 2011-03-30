#include "spamm_kernel.h"
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

/** Available stream kernels (their names). */
const char* spamm_stream_kernel_name [SPAMM_NUMBER_KERNELS] = {
  "kernel_experimental",
  "kernel_standard_SSE",
  "kernel_standard_SSE4_1",
  "kernel_Z_curve_SSE",
  "kernel_Z_curve_SSE4_1",
  "kernel_hierarchical_SSE",
  "kernel_hierarchical_SSE4_1",
  "kernel_standard_no_checks_SSE",
  "kernel_standard_no_checks_SSE4_1"
};

/** Get the name of the ith kernel.
 *
 * @param i The index of the kernel
 *
 * @return A pointer to a string that contains the name of the kernel.
 */
const char *
spamm_kernel_get_name (const unsigned int i)
{
  if (i >= SPAMM_NUMBER_KERNELS)
  {
    printf("illegal kernel index\n");
    exit(1);
  }

  return spamm_stream_kernel_name[i];
}

/** Return the kernel from a name.
 *
 * @param name The name of the kernel.
 *
 * @return The kernel.
 */
enum spamm_kernel_t
spamm_kernel_get_kernel (const char* name)
{
  enum spamm_kernel_t kernel;

  if (strcasecmp(name, "kernel_standard_SSE") == 0)
  {
    kernel = kernel_standard_SSE;
  }

  else if (strcasecmp(name, "kernel_standard_SSE4_1") == 0)
  {
    kernel = kernel_standard_SSE4_1;
  }

  else if (strcasecmp(name, "kernel_standard_no_checks_SSE") == 0)
  {
    kernel = kernel_standard_no_checks_SSE;
  }

  else if (strcasecmp(name, "kernel_standard_no_checks_SSE4_1") == 0)
  {
    kernel = kernel_standard_no_checks_SSE4_1;
  }

  else if (strcasecmp(name, "kernel_experimental") == 0)
  {
    kernel = kernel_experimental;
  }

  else if (strcasecmp(name, "kernel_Z_curve_SSE") == 0)
  {
    kernel = kernel_Z_curve_SSE;
  }

  else if (strcasecmp(name, "kernel_Z_curve_SSE4_1") == 0)
  {
    kernel = kernel_Z_curve_SSE4_1;
  }

  else if (strcasecmp(name, "kernel_hierarchical_SSE") == 0)
  {
    kernel = kernel_hierarchical_SSE;
  }

  else if (strcasecmp(name, "kernel_hierarchical_SSE4_1") == 0)
  {
    kernel = kernel_hierarchical_SSE4_1;
  }

  else
  {
    printf("unknown kernel: %s\n", name);
    exit(1);
  }

  return kernel;
}
