#ifndef __LOGGER_H
#define __LOGGER_H

#include "config.h"

#include "index.h"

#include <bitset>
#include <string>
#include <string.h>
#include <sstream>
#include <stdarg.h>
#include <stdio.h>

#ifdef DEBUG_OUTPUT
#define DEBUG(format, ...) logger(__FILE__, __LINE__, __func__, "", format, ##__VA_ARGS__)
#else
#define DEBUG(format, ...)
#endif

#define INFO(format, ...) logger(__FILE__, __LINE__, __func__, "INFO", format, ##__VA_ARGS__)
#define ABORT(format, ...) logger(__FILE__, __LINE__, __func__, "ERROR", format, ##__VA_ARGS__); CkExit()

#define FORMAT_LENGTH 2000
#define STRING_LENGTH 5000

inline
void logger (const char *const filename,
    const int linenumber,
    const char *const function_name,
    const char *const tag,
    const char *const format, ...)
{
  va_list ap;
  char new_format[FORMAT_LENGTH];
  char output_string[STRING_LENGTH];

  if(strlen(tag) > 0)
  {
    snprintf(new_format, FORMAT_LENGTH, "[%s:%d (%s) PE:%d %s] %s", filename,
        linenumber, function_name, CkMyPe(), tag, format);
  }

  else
  {
    snprintf(new_format, FORMAT_LENGTH, "[%s:%d (%s) PE:%d] %s", filename,
        linenumber, function_name, CkMyPe(), format);
  }

  va_start(ap, format);
  vsnprintf(output_string, STRING_LENGTH, new_format, ap);
  va_end(ap);
  CkPrintf(output_string);
}

/** Print a dense matrix.
 *
 * @param N The matrix size.
 * @param A The matrix.
 */
inline
void printDense (int N, double *A)
{
  std::ostringstream o;
  o.setf(std::ios::scientific);

  if(N <= 32)
  {
    for(int i = 0; i < N; i++) {
      for(int j = 0; j < N; j++)
      {
        o << " " << A[BLOCK_INDEX(i, j, 0, 0, N)];
      }
      o << std::endl;
    }
    CkPrintf(o.str().c_str());
  }

  else
  {
    INFO("matrix size too large for printing\n");
  }
}

inline
std::string toBinary (unsigned int i)
{
  std::string bitString = std::bitset<8*sizeof(unsigned int)>(i).to_string();
  while(bitString[0] == '0')
  {
    bitString.erase(0, 1);
  }
  return bitString;
}

#endif
