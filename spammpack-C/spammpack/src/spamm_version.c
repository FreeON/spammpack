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
  unsigned int length = strlen(PACKAGE_VERSION)+2+strlen(COMMIT_TAG)+1+1;
  char *version = (char*) malloc(sizeof(char)*length);

  snprintf(version, length, "%s (%s)", PACKAGE_VERSION, COMMIT_TAG);
  return version;
}
