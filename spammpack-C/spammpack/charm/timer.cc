#include <time.h>

#include "timer.h"

timer::timer ()
{
  timerRunning = false;
}

void timer::start ()
{
  clock_gettime(CLOCK_REALTIME, &startTime);
}

double timer::stop ()
{
  clock_gettime(CLOCK_REALTIME, &stopTime);
  timerRunning = false;

  return stopTime.tv_sec-startTime.tv_sec
    +(stopTime.tv_nsec-startTime.tv_nsec)/1.0e9;
}
