/** @file
 *
 * The implementation of the logging functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "logger.h"
#include "index.h"

#include <charm++.h>
#include <bitset>
#include <string>
#include <string.h>
#include <sstream>
#include <stdarg.h>
#include <stdio.h>


/** The main logging function. It is usually more convenient to use the macros
 * DEBUG, INFO, and ABORT instead.
 *
 * @param filename The name of the file where the message occurred.
 * @param linenumber The line number in that file.
 * @param function_name The name of the function.
 * @param tag A tag.
 * @param format The format string. See printf() for details.
 */
void logger (const char *const filename,
    const int linenumber,
    const char *const function_name,
    const char *const tag,
    const char *const format, ...)
{
  const int format_length = 2000;
  const int string_length = 5000;

  va_list ap;
  char new_format[format_length];
  char output_string[string_length];

  if(strlen(tag) > 0)
  {
    snprintf(new_format, format_length, "[%s:%d (%s) PE:%d %s] %s", filename,
        linenumber, function_name, CkMyPe(), tag, format);
  }

  else
  {
    snprintf(new_format, format_length, "[%s:%d (%s) PE:%d] %s", filename,
        linenumber, function_name, CkMyPe(), format);
  }

  va_start(ap, format);
  vsnprintf(output_string, string_length, new_format, ap);
  va_end(ap);
  CkPrintf(output_string);
}

/** Print a dense matrix.
 *
 * @param N The matrix size.
 * @param A The matrix.
 * @param format The format string for a message that is printed in front of
 * the matrix. See printf() for details.
 */
void printDense (int N, double *A, const char *const format, ...)
{
  std::ostringstream o;
  o.setf(std::ios::scientific);

  va_list ap;
  const int message_length = 2000;
  char message[message_length];

  va_start(ap, format);
  vsnprintf(message, message_length, format, ap);
  va_end(ap);

  o << message << std::endl;
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

/** Convert an integer to a binary as a string.
 *
 * @param i The integer.
 *
 * @return The string represenation.
 */
std::string toBinary (unsigned int i)
{
  std::string bitString = std::bitset<8*sizeof(unsigned int)>(i).to_string();
  while(bitString[0] == '0')
  {
    bitString.erase(0, 1);
  }
  return bitString;
}

/** Convert some number of bytes into a human readable form. For instance, n =
 * 1024 would result in "1 kiB".
 *
 * @param n The number of bytes to convert.
 *
 * @return The string representation of the human readable form.
 */
std::string humanReadableSize (unsigned long n)
{
  std::ostringstream o;
  double size = 0;
  const int numberUnits = 5;
  std::string unit[numberUnits] = { "B", "kiB", "MiB", "GiB", "TiB" };

  for(int i = 0; i < numberUnits; i++)
  {
    size /= 1024.;

    if(n/1024 > 0)
    {
      size += (n%1024);
      n /= 1024;
    }

    else
    {
      size += n;
      o << size << " " << unit[i];
      return o.str();
    }
  }
  ABORT("can not convert this number\n");
}
