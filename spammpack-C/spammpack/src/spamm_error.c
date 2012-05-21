#include "spamm_error.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void
spamm_error_fatal (const char *const filename, const int line, const char *const format, ...)
{
  va_list va;
  int format_length;
  char *new_format;

  format_length = 4+strlen(filename)+10+strlen(format);
  new_format = calloc(format_length, sizeof(char));
  snprintf(new_format, format_length, "[%s:%i] ", filename, line);
  strncat(new_format, format, format_length);

  /* Print error and exit. */
  va_start(va, format);
  vprintf(new_format, va);
  va_end(va);
  free(new_format);
  exit(1);
}
