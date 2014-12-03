/** @file */

#ifdef _OPENMP
#include <omp.h>
#endif

void
spamm_omp_init ()
{
#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_nested(1);
#endif
}
