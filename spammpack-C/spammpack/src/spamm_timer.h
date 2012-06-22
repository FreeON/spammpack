/** @file */

#ifndef __SPAMM_TIMER_H
#define __SPAMM_TIMER_H

#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

/** The timer type.
 */
enum spamm_timer_type_t
{
  /** Unitialized timer. */
  uninitialized,

  /** An empty timer. */
  empty,

  /** The walltime passed. */
  walltime,

  /** Total Instructions. */
  papi_total_instructions,

  /** Total cycles. */
  papi_total_cycles,

  /** Number floating point instructions. */
  papi_flop,

  /** Number single precision vector operations only. */
  papi_vec_sp,

  /** Number of L1 data-cache misses. */
  papi_l1_dcm,

  /** Number of L1 instruction-cache misses. */
  papi_l1_icm,

  /** Number of L2 data-cache misses. */
  papi_l2_dcm,

  /** Number of L2 instruction-cache misses. */
  papi_l2_icm
};

struct spamm_timer_t *
spamm_timer_new ();

void
spamm_timer_add_event (int event, struct spamm_timer_t *timer);

void
spamm_timer_delete (struct spamm_timer_t **timer);

void
spamm_timer_start (struct spamm_timer_t *timer);

void
spamm_timer_stop (struct spamm_timer_t *timer);

void
spamm_timer_get (short *length, unsigned long long **values, const struct spamm_timer_t *timer);

char *
spamm_timer_get_string (const struct spamm_timer_t *timer);

void
spamm_timer_info (const struct spamm_timer_t *timer, char *infostring, const int maxlength);

enum spamm_timer_type_t
spamm_timer_get_timer_type (const char *name);

void
spamm_timer_get_native_events ();

#endif
