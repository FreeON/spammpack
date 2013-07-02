#include <stdarg.h>
#include <stdio.h>
#include <string>
#include <sstream>

#include "Utilities.h"

Logging::Logging (const int loglevel)
{
  this->loglevel = loglevel;
}

void Logging::log (const int level,
    const char *const filename,
    const int line,
    const char *const format, ...)
{
  va_list va;
  std::string levelname;

  if(level >= Logging::DEBUG)
  {
    switch(level)
    {
      case Logging::DEBUG:
        levelname = "DEBUG";
        break;

      case Logging::INFO:
        levelname = "INFO";
        break;

      case Logging::ERROR:
        levelname = "ERROR";
        break;

      default:
        levelname = "UNKNOWN";
        break;
    };

    std::ostringstream format_stream;
    format_stream << "[" << filename << ":" << line << " - " << levelname << "] " << format;
    const char *new_format = format_stream.str().c_str();

    va_start(va, format);
    vprintf(new_format, va);
    va_end(va);
  }
}

void Logging::printBlock (const int level,
    const int chunksize,
    const std::string name,
    const float *const A)
{
  std::ostringstream debugString;
  debugString.precision(3);
  debugString.setf(std::ios::fixed);
  debugString << name << std::endl;
  for(int i = 0; i < chunksize; i++) {
    for(int j = 0; j < chunksize; j++)
    {
      debugString << " " << A[i+j*chunksize];
    }
    debugString << std::endl;
  }
  LOG_DEBUG(debugString.str().c_str());
}
