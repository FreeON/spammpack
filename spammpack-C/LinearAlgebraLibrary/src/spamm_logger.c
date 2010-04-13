#include "spamm.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void
spamm_log (const char *format, const char *filename, const int linenumber, ...)
{
  char *new_format;
  va_list ap;

  va_start(ap, linenumber);

  /* Fix up the format, i.e. insert file and line number. We conservatively
   * estimate that the file name and the line number will not exceed 2000
   * characters. */
  new_format = (char*) malloc(sizeof(char)*(strlen(format)+2000));
  sprintf(new_format, "[%s:%i] %s", filename, linenumber, format);

  vprintf(new_format, ap);

  free(new_format);
  va_end(ap);
}
