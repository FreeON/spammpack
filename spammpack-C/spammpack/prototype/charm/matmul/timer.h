#ifndef __TIMER_H
#define __TIMER_H

#include <string>
#include <time.h>

class Timer
{
  private:

    std::string message;
    struct timespec startTime;
    struct timespec endTime;

  public:

    Timer (std::string message);
    void stop ();
    const char * to_str ();
};

#endif
