#include "config.h"
#include "spamm.h"
#include <string.h>
#include <stdlib.h>

/** The version of the SpAMM library.
 *
 * @return The version.
 */
char *
spamm_version ()
{
  char *version = (char*) malloc(sizeof(char)*strlen(PACKAGE_VERSION));

  return version;
}
