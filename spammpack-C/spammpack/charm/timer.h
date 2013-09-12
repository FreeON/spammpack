/** @file
 *
 * The header file for the Timer class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __TIMER_H
#define __TIMER_H

#include <string>
#include <time.h>
#include <stdarg.h>

/** A timer object. */
class Timer
{
  private:

    /** A tag identifying this timer object. */
    std::string message;

    /** The start time. */
    struct timespec startTime;

    /** The end time. */
    struct timespec endTime;

    /** The string buffer for to_str(). */
    char *string_buffer;

  public:

    Timer (const char *format, ...);
    void stop ();
    const char * to_str ();
};

#endif
