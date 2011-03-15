#include "spamm_config.h"
#include "spamm_timer.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_PAPI
#include <papi.h>
#endif

/** The maximum length of the info output string. This is used only
 * internally. */
#define INFO_MAX_LENGTH 2000

/** @private A timer object.
 */
struct spamm_timer_t
{
  /** Indicate whether the timer is already running. */
  short timer_running;

  /** The timer type. */
  enum spamm_timer_type_t type;

#ifdef HAVE_PAPI
  /** The PAPI eventset. */
  int eventset;

  /** An array that holds the measured PAPI event values. */
  long_long *event_values;
#endif

  /** The time the timer was started. */
  struct rusage start_time;

  /** The time the timer was stopped. */
  struct rusage end_time;
};

#ifdef HAVE_PAPI
/** Handle PAPI related errors.
 *
 * @param error_code The error code returned by a PAPI call.
 * @param message The message to prepend to the error string.
 */
void
spamm_timer_handle_PAPI_error (const int error_code, const char *message)
{
  printf("[timer, PAPI error] %s: %s\n", message, PAPI_strerror(error_code));
  exit(1);
}
#endif

/** Return a new timer object.
 *
 * @param type The timer type.
 *
 * @return The newly allocated and initialized timer object.
 */
struct spamm_timer_t *
spamm_timer_new (const enum spamm_timer_type_t type)
{
#ifdef HAVE_PAPI
  int papi_result;
#endif
  struct spamm_timer_t *timer = (struct spamm_timer_t*) malloc(sizeof(struct spamm_timer_t));

  timer->timer_running = 0;
  timer->type = type;

#ifdef HAVE_PAPI
  timer->eventset = PAPI_NULL;

  if ((papi_result = PAPI_is_initialized()) != PAPI_LOW_LEVEL_INITED)
  {
    papi_result = PAPI_library_init(PAPI_VER_CURRENT);
    if (papi_result != PAPI_VER_CURRENT && papi_result > 0)
    {
      printf("PAPI library version mismatch\n");
      exit(1);
    }

    if (papi_result < 0)
    {
      spamm_timer_handle_PAPI_error(papi_result, "timer new, PAPI_library_init()");
    }

    if ((papi_result = PAPI_is_initialized()) != PAPI_LOW_LEVEL_INITED)
    {
      spamm_timer_handle_PAPI_error(papi_result, "timer new, PAPI_is_initialized()");
    }
  }
#endif

  switch (timer->type)
  {
    case walltime:
      break;

#ifdef HAVE_PAPI
    case papi_total_instructions:
    case papi_total_cycles:
    case papi_flop:
    case papi_vec_sp:
      if ((papi_result = PAPI_create_eventset(&timer->eventset)) != PAPI_OK)
      {
        spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_create_eventset()");
      }

      switch(timer->type)
      {
        case papi_total_instructions:
          if ((papi_result = PAPI_add_event(timer->eventset, PAPI_TOT_INS)) != PAPI_OK)
          {
            spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_add_event()");
          }

          timer->event_values = (long_long*) malloc(sizeof(long_long)*1);
          break;

        case papi_total_cycles:
          if ((papi_result = PAPI_add_event(timer->eventset, PAPI_TOT_CYC)) != PAPI_OK)
          {
            spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_add_event()");
          }

          timer->event_values = (long_long*) malloc(sizeof(long_long)*1);
          break;

        case papi_flop:
          if ((papi_result = PAPI_add_event(timer->eventset, PAPI_FP_OPS)) != PAPI_OK)
          {
            spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_add_event()");
          }

          if ((papi_result = PAPI_add_event(timer->eventset, PAPI_TOT_CYC)) != PAPI_OK)
          {
            spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_add_event()");
          }

          timer->event_values = (long_long*) malloc(sizeof(long_long)*2);
          break;

        case papi_vec_sp:
          if ((papi_result = PAPI_add_event(timer->eventset, PAPI_VEC_SP)) != PAPI_OK)
          {
            spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_add_event()");
          }

          timer->event_values = (long_long*) malloc(sizeof(long_long)*1);
          break;
      }
      break;
#endif

    default:
      printf("[timer new] unknown timer type\n");
      exit(1);
      break;
  }

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
#ifdef HAVE_PAPI
  int papi_result;
#endif

  switch ((*timer)->type)
  {
    case walltime:
      break;

#ifdef HAVE_PAPI
    case papi_total_instructions:
    case papi_total_cycles:
    case papi_flop:
    case papi_vec_sp:
      if ((papi_result = PAPI_cleanup_eventset((*timer)->eventset)) != PAPI_OK)
      {
        spamm_timer_handle_PAPI_error(papi_result, "delete timer, PAPI_cleanup_eventset()");
      }

      if ((papi_result = PAPI_destroy_eventset(&(*timer)->eventset)) != PAPI_OK)
      {
        spamm_timer_handle_PAPI_error(papi_result, "delete timer, PAPI_destroy_eventset()");
      }

      free((*timer)->event_values);
      break;
#endif

    default:
      printf("[timer delete] unknown timer type\n");
      exit(1);
      break;
  }

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
#ifdef HAVE_PAPI
  int papi_result;
#endif

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

#ifdef HAVE_PAPI
    case papi_total_instructions:
    case papi_total_cycles:
    case papi_flop:
    case papi_vec_sp:
      if ((papi_result = PAPI_start(timer->eventset)) != PAPI_OK)
      {
        spamm_timer_handle_PAPI_error(papi_result, "start timer, PAPI_start()");
      }
      break;
#endif

    default:
      printf("[timer start] unknown timer type\n");
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
#ifdef HAVE_PAPI
  int papi_result;
#endif

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

#ifdef HAVE_PAPI
    case papi_total_instructions:
    case papi_total_cycles:
    case papi_flop:
    case papi_vec_sp:
      if ((papi_result = PAPI_stop(timer->eventset, timer->event_values)) != PAPI_OK)
      {
        spamm_timer_handle_PAPI_error(papi_result, "stop timer, PAPI_stop()");
      }
      break;
#endif

    default:
      printf("[timer stop] unknown timer type\n");
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
unsigned long long
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

#ifdef HAVE_PAPI
    case papi_total_instructions:
    case papi_total_cycles:
    case papi_flop:
    case papi_vec_sp:
      return timer->event_values[0];
      break;
#endif

    default:
      printf("[timer get] unknown timer type\n");
      exit(1);
      break;
  }
}

/** Get the floprate. This is only possible for timer type papi_flop.
 *
 * @param timer Tht timer.
 *
 * @return The flop rate in Mflop/s.
 */
float
spamm_timer_get_floprate (const struct spamm_timer_t *timer)
{
#ifdef HAVE_PAPI
  int papi_result;
#endif

  if (timer->timer_running != 0)
  {
    printf("[timer get] this timer is running right now\n");
    exit(1);
  }

  switch (timer->type)
  {
#ifdef HAVE_PAPI
    case papi_flop:
      if ((papi_result = PAPI_get_opt(PAPI_CLOCKRATE, NULL)) <= 0)
      {
        spamm_timer_handle_PAPI_error(papi_result, "[timer get] failed PAPI_get_opt()");
      }
      return (float) timer->event_values[0] / (float) timer->event_values[1] * (float) papi_result;
      break;
#endif
    default:
      printf("[timer get] unknown or illegal timer type\n");
      exit(1);
      break;
  }
}

/** Print out some information on the timer.
 *
 * @param timer The timer to print information on.
 * @param infostring A pointer to the output string. This has to be allocated
 * by the caller.
 * @param maxlength The length of the output string.
 */
void
spamm_timer_info (const struct spamm_timer_t *timer, char *infostring,
    const int maxlength)
{
  int i;
  char string[INFO_MAX_LENGTH];

#ifdef HAVE_PAPI
  int papi_result;
#endif

  assert(timer != NULL);
  assert(infostring != NULL);
  assert(maxlength >= 1);

  /* Terminate string. */
  string[0] = '\0';

  switch (timer->type)
  {
    case walltime:
      sprintf(string, "walltime");
      break;

#ifdef HAVE_PAPI
    case papi_total_instructions:
      sprintf(string, "PAPI Total instructions");
      break;

    case papi_total_cycles:
      if ((papi_result = PAPI_get_opt(PAPI_CLOCKRATE, NULL)) <= 0)
      {
        spamm_timer_handle_PAPI_error(papi_result, "timer info, PAPI_get_opt()");
      }
      sprintf(string, "PAPI Total cycles, clockrate = %d MHz", papi_result);
      break;

    case papi_flop:
      if ((papi_result = PAPI_get_opt(PAPI_CLOCKRATE, NULL)) <= 0)
      {
        spamm_timer_handle_PAPI_error(papi_result, "timer info, PAPI_get_opt()");
      }
      sprintf(string, "PAPI Floating point operations, clockrate = %d MHz", papi_result);
      break;

    case papi_vec_sp:
      sprintf(string, "PAPI Single precision vector/SIMD instructions");
      break;
#endif

    default:
      printf("[timer info] unknown timer type\n");
      exit(1);
      break;
  }

  /* Copy string to output. */
  for (i = 0; i < maxlength && i < INFO_MAX_LENGTH; i++)
  {
    infostring[i] = string[i];
    if (string[i] == '\0')
    {
      break;
    }
  }
}
