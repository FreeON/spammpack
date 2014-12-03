/** @file */

#define CONCAT2(a, b) a ## _ ## b
#define SPAMM_FUNCTION(name, type) CONCAT2(name, type)

#define FUNCNAME sgemm
#define FUNCTYPE float
#include "spamm_blas_source.c"
#undef FUNCTYPE
#undef FUNCNAME

#define FUNCNAME dgemm
#define FUNCTYPE double
#include "spamm_blas_source.c"
#undef FUNCTYPE
#undef FUNCNAME
