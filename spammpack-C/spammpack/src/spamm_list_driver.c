/** This code generates all types of spamm_list.h
 */

/** This is kind of a hack. In order to expand TYPENAME it seems that we need
 * to pass it twice as an argument.
 */
#define CONCAT3_REAL(a, b, c) a ## _ ## b ## _ ## c
#define CONCAT3(a, b, c) CONCAT3_REAL(a, b, c)

#define SYMBOL(symbolname) CONCAT3(spamm, TYPENAME, symbolname)

#define TYPE unsigned int
#define TYPENAME unsigned_int
#include "spamm_list_source.c"
#undef TYPENAME
#undef TYPE

#define TYPE float
#define TYPENAME float
#include "spamm_list_source.c"
#undef TYPENAME
#undef TYPE
