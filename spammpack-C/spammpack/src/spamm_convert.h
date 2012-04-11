#ifndef __SPAMM_CONVERT_H
#define __SPAMM_CONVERT_H

#include "spamm_types_private.h"

struct spamm_t *
spamm_convert_dense_to_spamm (const unsigned int M, const unsigned int N,
    const enum spamm_layout_t dense_type, float *A_dense,
    const enum spamm_layout_t spamm_layout);

#endif
