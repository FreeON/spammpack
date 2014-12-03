#ifndef __UTILITIES_H
#define __UTILITIES_H

#define LOG_DEBUG(format, ...) Logging::log(Logging::DEBUG, __FILE__, __LINE__, format, ##__VA_ARGS__)
#define LOG_INFO(format, ...)  Logging::log(Logging::INFO,  __FILE__, __LINE__, format, ##__VA_ARGS__)
#define LOG_ERROR(format, ...) Logging::log(Logging::ERROR, __FILE__, __LINE__, format, ##__VA_ARGS__)

#include <string>

class Logging
{
  private:

    int loglevel;

  public:

    static const int DEBUG = 0;
    static const int INFO  = 1;
    static const int ERROR = 2;

    Logging (const int loglevel);

    static void log (const int level, const char *const filename,
        const int line, const char *const format, ...);
    static void printBlock (const int level, const int chunksize,
        const std::string name, const float *const A);
};

class Counter
{
  public:

    static void reset ();
    static void increment ();
    static int get ();
};

#endif
