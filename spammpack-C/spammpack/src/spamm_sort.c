/** @file */

#define SPAMM_FUNCTION(base,type) base##type

#define SPAMM_SORT_TYPE unsigned int
#define SPAMM_FUNC_TYPE unsigned_int
#include "spamm_sort_source.c"
#undef SPAMM_SORT_TYPE
#undef SPAMM_FUNC_TYPE
