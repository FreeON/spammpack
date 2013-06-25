#ifndef __TIMER_H
#define __TIMER_H

class timer
{
  private:

    bool timerRunning;

    struct timespec startTime;
    struct timespec stopTime;

  public:

    timer ();
    void start ();
    double stop ();
};

#endif
