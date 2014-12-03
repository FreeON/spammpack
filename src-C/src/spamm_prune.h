/** @file */

#ifndef __SPAMM_PRUNE_H
#define __SPAMM_PRUNE_H

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
spamm_prune (struct spamm_matrix_t *const A);

__END_DECLARATIONS

#endif
