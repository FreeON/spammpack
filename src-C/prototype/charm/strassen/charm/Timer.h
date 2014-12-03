#ifndef __TIMER_H
#define __TIMER_H

#include <string>

class Timer
{
  private:

    std::string name;
    bool isRunning;
    struct timespec startTime;
    struct timespec endTime;

  public:

    Timer (std::string name);
    void start ();
    void stop ();
};

#endif
