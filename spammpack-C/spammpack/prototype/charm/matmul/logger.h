#ifndef __LOGGER_H
#define __LOGGER_H

#include "config.h"

#include <stdarg.h>

#ifdef DEBUG_OUTPUT
#define DEBUG(format, ...) logger(__FILE__, __LINE__, "", format, ##__VA_ARGS__)
#else
#define DEBUG(format, ...)
#endif

#define ABORT(format, ...) logger(__FILE__, __LINE__, "ERROR", format, ##__VA_ARGS__); CkExit()

#define FORMAT_LENGTH 2000
#define STRING_LENGTH 5000

inline
void logger (const char *const filename, const int linenumber,
    const char *const tag,
    const char *const format, ...)
{
  va_list ap;
  char new_format[FORMAT_LENGTH];
  char output_string[STRING_LENGTH];

  if(strlen(tag) > 0)
  {
    snprintf(new_format, FORMAT_LENGTH, "[%s:%d PE:%d %s] %s", filename, linenumber, CkMyPe(), tag, format);
  }

  else
  {
    snprintf(new_format, FORMAT_LENGTH, "[%s:%d PE:%d] %s", filename, linenumber, CkMyPe(), format);
  }

  va_start(ap, format);
  vsnprintf(output_string, STRING_LENGTH, new_format, ap);
  va_end(ap);
  CkPrintf(output_string);
}

#endif
