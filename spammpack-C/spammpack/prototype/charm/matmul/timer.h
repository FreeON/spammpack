#ifndef __TIMER_H
#define __TIMER_H

#include <string>
#include <time.h>
#include <stdarg.h>

class Timer
{
  private:

    std::string message;
    struct timespec startTime;
    struct timespec endTime;

  public:

    Timer (const char *format, ...);
    void stop ();
    const char * to_str ();
};

#endif
