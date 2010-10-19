#include "config.h"
#include "spamm.h"
#include <string.h>
#include <stdlib.h>

/** The version of the SpAMM library.
 *
 * @return A string is returned that holds the version of the library. This
 * string must be freed by the caller.
 */
char *
spamm_version ()
{
  char *version = (char*) malloc(sizeof(char)*(strlen(PACKAGE_VERSION)+1));

  strncpy(version, PACKAGE_VERSION, strlen(PACKAGE_VERSION)+1);
  return version;
}
