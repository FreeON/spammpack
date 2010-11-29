#ifndef __SPAMM_TIMER_H
#define __SPAMM_TIMER_H

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

/** A timer object.
 */
struct spamm_timer_t
{
  short timer_running;
  struct rusage start_time;
  struct rusage end_time;
};

struct spamm_timer_t *
spamm_timer_new ();

void
spamm_timer_delete (struct spamm_timer_t **timer);

void
spamm_timer_start (struct spamm_timer_t *timer);

void
spamm_timer_stop (struct spamm_timer_t *timer);

float
spamm_timer_get_seconds (const struct spamm_timer_t *timer);

#endif
