#include "spamm_timer.h"
#include <stdlib.h>
#include <stdio.h>

/** Return a new timer object.
 *
 * @return The newly allocated and initialized timer object.
 */
struct spamm_timer_t *
spamm_timer_new ()
{
  struct spamm_timer_t *timer = (struct spamm_timer_t*) malloc(sizeof(struct spamm_timer_t));

  timer->timer_running = 0;

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
  if (getrusage(RUSAGE_SELF, &timer->start_time) != 0)
  {
    printf("[start timer] error getting time\n");
    exit(1);
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
  if (getrusage(RUSAGE_SELF, &timer->end_time) != 0)
  {
    printf("[stop timer] error getting time\n");
    exit(1);
  }
}

/** Get the number of seconds that have passed between starting and stopping
 * the timer. This time is using the total time, i.e. the sum of user and
 * system time.
 *
 * @param timer The timer.
 */
float
spamm_timer_get_seconds (const struct spamm_timer_t *timer)
{
  unsigned int seconds, useconds;

  if (timer->timer_running != 0)
  {
    printf("[get seconds] this timer is running right now\n");
    exit(1);
  }

  seconds = (timer->end_time).ru_utime.tv_sec
    +(timer->end_time).ru_stime.tv_sec
    -((timer->start_time).ru_utime.tv_sec
        +(timer->start_time).ru_stime.tv_sec);

  useconds = (timer->end_time).ru_utime.tv_usec
    +(timer->end_time).ru_stime.tv_usec
    -((timer->start_time).ru_utime.tv_usec
        +(timer->start_time).ru_stime.tv_usec);

  return (double) seconds + (double) useconds/1.0e6;
}
