/** @file */

#ifndef __SPAMM_SORT_H
#define __SPAMM_SORT_H

#include "spamm_types.h"

#ifdef __cplusplus
#define __BEGIN_DECLARATIONS extern "C" {
#define __END_DECLARATIONS }
#else
#define __BEGIN_DECLARATIONS
#define __END_DECLARATIONS
#endif

__BEGIN_DECLARATIONS

void
spamm_sort_masked (const unsigned int length,
    unsigned int *list,
    const unsigned int mask);

void
spamm_sort_norm (const unsigned int length,
    unsigned int *list,
    spamm_norm_t *norm);

__END_DECLARATIONS

#endif
