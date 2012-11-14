/** @file */

#define CONCAT2(a, b) a ## _ ## b
#define SPAMM_FUNCTION(name, type) CONCAT2(name, type)

#define SPAMM_SORT_TYPE unsigned int
#define SPAMM_FUNC_TYPE unsigned_int
#define SPAMM_SORT_MASKED
#define SPAMM_SORT_COMPARE(a, b, mask) ((a) & mask < (b) & mask ? -1 : 1)
#include "spamm_sort_source.c"
#undef SPAMM_SORT_TYPE
#undef SPAMM_FUNC_TYPE
#undef SPAMM_SORT_MASKED
#undef SPAMM_SORT_COMPARE

#define SPAMM_SORT_TYPE float
#define SPAMM_FUNC_TYPE float
#include "spamm_sort_source.c"
#undef SPAMM_SORT_TYPE
#undef SPAMM_FUNC_TYPE
