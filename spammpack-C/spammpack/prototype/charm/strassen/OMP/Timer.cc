#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "strassenOMP.h"

Timer::Timer (std::string name)
{
  isRunning = false;
  this->name = name;
}

void Timer::start ()
{
  if(isRunning)
  {
    printf("timer already running\n");
    exit(1);
  }

  if(clock_gettime(CLOCK_MONOTONIC, &startTime) < 0)
  {
    printf("error starting timer\n");
    exit(1);
  }

  isRunning = true;
  printf("starting %s... ", name.c_str()); fflush(stdout);
}

void Timer::stop ()
{
  if(!isRunning)
  {
    printf("timer is not running\n");
    exit(1);
  }

  if(clock_gettime(CLOCK_MONOTONIC, &endTime) < 0)
  {
    printf("error stopping timer\n");
    exit(1);
  }

  isRunning = false;
  printf("%f seconds\n", endTime.tv_sec+endTime.tv_nsec/1.0e9-(startTime.tv_sec+startTime.tv_nsec/1.0e9));
}
