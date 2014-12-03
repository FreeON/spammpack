/** @file */

#define CONCAT2(a, b) a ## _ ## b
#define SPAMM_FUNCTION(name, type) CONCAT2(name, type)

#define SPAMM_SORT_MASKED
#include "spamm_sort_source.c"
#undef SPAMM_SORT_MASKED

#define SPAMM_SORT_NORM
#include "spamm_sort_source.c"
#undef SPAMM_SORT_NORM
