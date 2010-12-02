/** @file */

#ifndef __SPAMM_TIMER_H
#define __SPAMM_TIMER_H

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

/** A timer object.
 */
struct spamm_timer_t
{
  /** Indicate whether the timer is already running. */
  short timer_running;

  /** The time the timer was started. */
  struct rusage start_time;

  /** The time the timer was stopped. */
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
