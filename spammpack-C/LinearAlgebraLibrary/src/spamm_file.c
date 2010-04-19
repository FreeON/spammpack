#include "spamm.h"

#if defined(HAVE_BZLIB_H) && defined(HAVE_LIBBZ2)
#include <bzlib.h>

void
spamm_bz2_open (const char *filename)
{
}

#else

void
spamm_bz2_open (const char *filename)
{
}

#endif
