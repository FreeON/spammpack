#include "timer.h"
#include <string>
#include <sstream>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

Timer::Timer (std::string message)
{
  this->message = message;
  if(clock_gettime(CLOCK_MONOTONIC_RAW, &startTime) < 0)
  {
    printf("can not start timer\n");
    exit(1);
  }
}

void Timer::stop ()
{
  if(clock_gettime(CLOCK_MONOTONIC_RAW, &endTime) < 0)
  {
    printf("can not stop timer\n");
    exit(1);
  }
}

const char * Timer::to_str ()
{
  std::ostringstream o;
  o.setf(std::ios::fixed);
  o << message << ": ";
  o << endTime.tv_sec+endTime.tv_nsec/1.0e9
    -(startTime.tv_sec+startTime.tv_nsec/1.0e9);
  o << " seconds" << std::endl;
  return o.str().c_str();
}
