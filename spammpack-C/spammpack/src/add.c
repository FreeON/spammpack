#include "typed_function.h"

#define SPAMM_TYPE float
#include "add_source.c"
#undef SPAMM_TYPE

#define SPAMM_TYPE double
#include "add_source.c"
#undef SPAMM_TYPE
