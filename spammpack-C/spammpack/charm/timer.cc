/** @file
 *
 * Implementation of the Timer class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"
#include "timer.h"
#include <string>
#include <string.h>
#include <sstream>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

/** The buffer size for string outputs. */
#define BUFFER_SIZE 2000

/** The constructor.
 *
 * @param format The format string. See printf() for details.
 */
Timer::Timer (const char *format, ...)
{
  char output_buffer[BUFFER_SIZE];
  va_list ap;

  va_start(ap, format);
  vsnprintf(output_buffer, BUFFER_SIZE, format, ap);
  va_end(ap);

  startTime.tv_sec = 0;
  startTime.tv_nsec = 0;

  endTime.tv_sec = 0;
  endTime.tv_nsec = 0;

  message = std::string(output_buffer);
  string_buffer = NULL;

  timerRunning = false;

  //printf("new timer %s\n", message.c_str());
}

/** The destructor. */
Timer::~Timer (void)
{
  if(string_buffer != NULL)
  {
    free(string_buffer);
  }

  //printf("deleting timer %s\n", message.c_str());
}

/** Reset the timer and start. */
void Timer::start (void)
{
  if(timerRunning)
  {
    printf("[%s:%d] timer is already running\n", __FILE__, __LINE__);
    exit(1);
  }

  if(clock_gettime(CLOCKTYPE, &startTime) < 0)
  {
    printf("[%s;%d] can not start timer\n", __FILE__, __LINE__);
    exit(1);
  }

  timerRunning = true;
}

/** Stop the timer. */
void Timer::stop (void)
{
  if(!timerRunning)
  {
    printf("[%s:%d] timer is not running\n", __FILE__, __LINE__);
    exit(1);
  }

  if(clock_gettime(CLOCKTYPE, &endTime) < 0)
  {
    printf("[%s;%d] can not stop timer\n", __FILE__, __LINE__);
    exit(1);
  }

  timerRunning = false;
}

/** Convert the timer to a string.
 *
 * @return The timer as a string. The string should not be free()'ed.
 */
const char * Timer::to_str (void)
{
  std::ostringstream o;

  if(timerRunning)
  {
    printf("[%s;%d] timer is running\n", __FILE__, __LINE__);
    exit(1);
  }

  o.setf(std::ios::fixed);
  o << message << ": ";
  o << endTime.tv_sec+endTime.tv_nsec/1.0e9
    -(startTime.tv_sec+startTime.tv_nsec/1.0e9);
  o << " seconds";
  if(string_buffer != NULL)
  {
    free(string_buffer);
  }
  string_buffer = strdup(o.str().c_str());
  return string_buffer;
}
