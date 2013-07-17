#ifndef __LOG_H
#define __LOG_H

#include "config.h"

#include <stdarg.h>

#ifdef DEBUG
#define LOG(format, ...) logger(__FILE__, __LINE__, "", format, ##__VA_ARGS__)
#else
#define LOG(format, ...)
#endif

#define ERROR(format, ...) logger(__FILE__, __LINE__, "ERROR", format, ##__VA_ARGS__)

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
    snprintf(new_format, FORMAT_LENGTH, "[%s:%d %s] %s", filename, linenumber, tag, format);
  }

  else
  {
    snprintf(new_format, FORMAT_LENGTH, "[%s:%d] %s", filename, linenumber, format);
  }

  va_start(ap, format);
  vsnprintf(output_string, STRING_LENGTH, new_format, ap);
  va_end(ap);
  CkPrintf(output_string);
}

#endif
