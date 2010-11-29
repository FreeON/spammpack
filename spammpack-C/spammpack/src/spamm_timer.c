#include "spamm_timer.h"
#include <stdlib.h>
#include <stdio.h>

struct spamm_timer_t *
spamm_timer_new ()
{
  struct spamm_timer_t *timer = (struct spamm_timer_t*) malloc(sizeof(struct spamm_timer_t));

  timer->timer_running = 0;

  return timer;
}

void
spamm_timer_delete (struct spamm_timer_t **timer)
{
  free(*timer);
  *timer = NULL;
}

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

float
spamm_timer_get_seconds (const struct spamm_timer_t *timer)
{
  if (timer->timer_running != 0)
  {
    printf("[print timer] this timer is running right now\n");
    exit(1);
  }

  return (timer->end_time).ru_utime.tv_sec+(timer->end_time).ru_utime.tv_usec/1.0e6
    -((timer->start_time).ru_utime.tv_sec+(timer->start_time).ru_utime.tv_usec/1.0e6);
}
