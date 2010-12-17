#include "config.h"
#include "spamm_timer.h"
#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_PAPI
#include <papi.h>
#endif

/** @private A timer object.
 */
struct spamm_timer_t
{
  /** Indicate whether the timer is already running. */
  short timer_running;

  /** The timer type. */
  enum spamm_timer_type_t type;

#ifdef HAVE_PAPI
#endif

  /** The time the timer was started. */
  struct rusage start_time;

  /** The time the timer was stopped. */
  struct rusage end_time;
};

/** Return a new timer object.
 *
 * @return The newly allocated and initialized timer object.
 */
struct spamm_timer_t *
spamm_timer_new (const enum spamm_timer_type_t type)
{
  struct spamm_timer_t *timer = (struct spamm_timer_t*) malloc(sizeof(struct spamm_timer_t));

  timer->timer_running = 0;
  timer->type = type;

  return timer;
}

/** Delete a timer object.
 *
 * @param timer The timer object to delete. The memory pointer to will be
 * freed and the pointer to the timer will be set to NULL.
 */
void
spamm_timer_delete (struct spamm_timer_t **timer)
{
  free(*timer);
  *timer = NULL;
}

/** Start the timer.
 *
 * @param timer The timer to start.
 */
void
spamm_timer_start (struct spamm_timer_t *timer)
{
  if (timer->timer_running != 0)
  {
    printf("[start timer] this timer is already running\n");
    exit(1);
  }

  timer->timer_running = 1;

  switch (timer->type)
  {
    case walltime:
      if (getrusage(RUSAGE_SELF, &timer->start_time) != 0)
      {
        printf("[start timer] error getting time\n");
        exit(1);
      }
      break;

    default:
      printf("[FIXME] unknown timer type\n");
      exit(1);
      break;
  }
}

/** Stop the timer.
 *
 * @param timer The timer to stop.
 */
void
spamm_timer_stop (struct spamm_timer_t *timer)
{
  if (timer->timer_running != 1)
  {
    printf("[stop timer] this timer is not running\n");
    exit(1);
  }

  timer->timer_running = 0;

  switch (timer->type)
  {
    case walltime:
      if (getrusage(RUSAGE_SELF, &timer->end_time) != 0)
      {
        printf("[stop timer] error getting time\n");
        exit(1);
      }
      break;

    default:
      printf("[FIXME] unknown timer type\n");
      exit(1);
      break;
  }
}

/** Get the time passed between starting and stopping the timer. This time is
 * is returned as an integer. The units of that integer depends on the timer
 * type.
 *
 * @param timer The timer.
 *
 * @return The time passed in the default units of the timer type.
 */
unsigned int
spamm_timer_get (const struct spamm_timer_t *timer)
{
  if (timer->timer_running != 0)
  {
    printf("[timer get] this timer is running right now\n");
    exit(1);
  }

  switch (timer->type)
  {
    case walltime:
      return 1000000*((timer->end_time).ru_utime.tv_sec
          +(timer->end_time).ru_stime.tv_sec
          -((timer->start_time).ru_utime.tv_sec
            +(timer->start_time).ru_stime.tv_sec))
        +(timer->end_time).ru_utime.tv_usec
        +(timer->end_time).ru_stime.tv_usec
        -((timer->start_time).ru_utime.tv_usec
            +(timer->start_time).ru_stime.tv_usec);
      break;

    default:
      printf("[FIXME] unknown timer type\n");
      exit(1);
      break;
  }
}
