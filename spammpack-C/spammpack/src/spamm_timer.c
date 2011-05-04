#include "config.h"
#include "spamm_timer.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

#ifdef HAVE_PAPI
#include <papi.h>
#endif

/** @private A timer object.
 */
struct spamm_timer_t
{
  /** Indicate whether the timer is already running. */
  short timer_running;

#ifdef HAVE_PAPI
  /** The PAPI eventset. */
  int eventset;

  /** The length of the values array. */
  short number_events;

  /** An array holding the PAPI event codes. */
  int *event;

  /** An array that holds the measured PAPI event values. */
  long_long *event_value;
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

/** Initialize the PAPI library.
 */
void
spamm_timer_PAPI_init ()
{
#ifdef HAVE_PAPI
  int papi_result;

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
#else
  printf("PAPI is not available\n");
  exit(1);
#endif
}

/** Return a new timer object.
 *
 * @return The newly allocated and initialized timer object.
 */
struct spamm_timer_t *
spamm_timer_new ()
{
#ifdef HAVE_PAPI
  int papi_result;
#endif

  struct spamm_timer_t *timer = calloc(1, sizeof(struct spamm_timer_t));

  timer->timer_running = 0;

#ifdef HAVE_PAPI
  timer->eventset = PAPI_NULL;
  spamm_timer_PAPI_init();
  if ((papi_result = PAPI_create_eventset(&timer->eventset)) != PAPI_OK)
  {
    spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_create_eventset()");
  }
#endif

  return timer;
}

/** Add an event to the timer.
 *
 * @param event The event to add.
 * @param timer The timer to add to.
 */
void
spamm_timer_add_event (int event, struct spamm_timer_t *timer)
{
#ifdef HAVE_PAPI
  short i;
  int papi_result;
  int *temp_event;
  char temp_string[100];

  if ((papi_result = PAPI_add_event(timer->eventset, event)) != PAPI_OK)
  {
    spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_add_event()");
  }

  temp_event = calloc(timer->number_events+1, sizeof(int));
  for (i = 0; i < timer->number_events; i++)
  {
    temp_event[i] = timer->event[i];
  }
  temp_event[timer->number_events] = event;

  if (timer->event != NULL)
  {
    free(timer->event);
  }
  timer->event = temp_event;
  timer->number_events += 1;
#else

  printf("[add event] adding walltime event\n");

#endif
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
  int i;
  int papi_result;

  if ((papi_result = PAPI_cleanup_eventset((*timer)->eventset)) != PAPI_OK)
  {
    spamm_timer_handle_PAPI_error(papi_result, "delete timer, PAPI_cleanup_eventset()");
  }

  if ((papi_result = PAPI_destroy_eventset(&(*timer)->eventset)) != PAPI_OK)
  {
    spamm_timer_handle_PAPI_error(papi_result, "delete timer, PAPI_destroy_eventset()");
  }
  free((*timer)->event);
  free((*timer)->event_value);
#endif

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

#ifdef HAVE_PAPI
  if (timer->number_events == 0)
  {
    printf("[start timer] this timer does not have any events\n");
    exit(1);
  }

  if ((papi_result = PAPI_start(timer->eventset)) != PAPI_OK)
  {
    spamm_timer_handle_PAPI_error(papi_result, "start timer, PAPI_start()");
  }
#else
  if (getrusage(RUSAGE_SELF, &timer->start_time) != 0)
  {
    printf("[start timer] error getting time\n");
    exit(1);
  }
#endif
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

#ifdef HAVE_PAPI
  if (timer->event_value != NULL)
  {
    free(timer->event_value);
  }
  timer->event_value = calloc(timer->number_events, sizeof(unsigned long long));
  if ((papi_result = PAPI_stop(timer->eventset, timer->event_value)) != PAPI_OK)
  {
    spamm_timer_handle_PAPI_error(papi_result, "stop timer, PAPI_stop()");
  }
#else
  if (getrusage(RUSAGE_SELF, &timer->end_time) != 0)
  {
    printf("[stop timer] error getting time\n");
    exit(1);
  }
#endif
}

/** Get the time passed between starting and stopping the timer. This time is
 * is returned as an integer. The units of that integer depends on the timer
 * type.
 *
 * @param length Will hold the length of the allocated array values on return.
 * @param timer The timer.
 *
 * @return The variables length will contain the number of elements in the
 * values array, and the values variable will point to an array with the
 * counter values. This array needs to be free()'ed by the caller when not
 * needed anymore.
 */
void
spamm_timer_get (short *length, unsigned long long **values, const struct spamm_timer_t *timer)
{
#ifdef HAVE_PAPI
  short i;
#endif

  assert(length != NULL);
  assert(values != NULL);

  if (timer->timer_running != 0)
  {
    printf("[timer get] this timer is running right now\n");
    exit(1);
  }

#ifdef HAVE_PAPI
  *length = timer->number_events;
  *values = calloc(*length, sizeof(unsigned long long));
  for (i = 0; i < *length; i++)
  {
    (*values)[i] = timer->event_value[i];
  }
#else
  *length = 1;
  *values = calloc(*length, sizeof(unsigned long long));
  (*values)[0] = 1000000*((timer->end_time).ru_utime.tv_sec
      +(timer->end_time).ru_stime.tv_sec
      -((timer->start_time).ru_utime.tv_sec
        +(timer->start_time).ru_stime.tv_sec))
    +(timer->end_time).ru_utime.tv_usec
    +(timer->end_time).ru_stime.tv_usec
    -((timer->start_time).ru_utime.tv_usec
        +(timer->start_time).ru_stime.tv_usec);
#endif
}

/** Return a string with all timer values in list format.
 *
 * timer The timer.
 *
 * @return A string with a list of timer values. The string needs to be
 * free()'ed by the caller.
 */
char *
spamm_timer_get_string (const struct spamm_timer_t *timer)
{
  const int maxlength = 2000;
  int i;
  char *result = calloc(maxlength, sizeof(char));

#ifdef HAVE_PAPI
  int papi_result;
  char *temp_string = calloc(maxlength, sizeof(char));
  char *temp_event_name = calloc(maxlength, sizeof(char));

  sprintf(result, "[");
  for (i = 0; i < timer->number_events; i++)
  {
    if ((papi_result = PAPI_event_code_to_name(timer->event[i], temp_event_name)) != PAPI_OK)
    {
      spamm_timer_handle_PAPI_error(papi_result, "get string, PAPI_event_code_to_name()");
    }

    sprintf(temp_string, " %lli 0x%x(%s)", timer->event_value[i], timer->event[i], temp_event_name);
    strncat(result, temp_string, maxlength-1);
    if (i < timer->number_events-1)
    {
      strncat(result, ",", maxlength-1);
    }
  }
  strncat(result, " ]", maxlength-1);
  free(temp_string);
#else
  sprintf(result, "%u (walltime)", 0);
#endif

  return result;
}

/** Print out some information on the timer.
 *
 * @param timer The timer to print information on.
 * @param infostring A pointer to the output string. This has to be allocated
 * by the caller.
 * @param maxlength The length of the output string.
 */
void
spamm_timer_info (const struct spamm_timer_t *timer, char *infostring, const int maxlength)
{
  int i;

#ifdef HAVE_PAPI
  int papi_result;
  char *temp_string;
  char *temp_event_name;
#endif

  assert(timer != NULL);
  assert(infostring != NULL);
  assert(maxlength >= 1);

#ifdef HAVE_PAPI
  temp_string = calloc(maxlength, sizeof(char));
  temp_event_name = calloc(maxlength, sizeof(char));

  sprintf(infostring, "[");
  for (i = 0; i < timer->number_events; i++)
  {
    if ((papi_result = PAPI_event_code_to_name(timer->event[i], temp_event_name)) != PAPI_OK)
    {
      spamm_timer_handle_PAPI_error(papi_result, "timer info, PAPI_event_code_to_name()");
    }

    sprintf(temp_string, " 0x%x(%s)", timer->event[i], temp_event_name);
    strncat(infostring, temp_string, maxlength-1);
    if (i < timer->number_events-1)
    {
      strncat(infostring, ",", maxlength-1);
    }
  }
  strncat(infostring, " ]", maxlength-1);
  free(temp_string);
#else
  sprintf(infostring, "walltime");
#endif
}

/** Get timer from a name.
 *
 * @name The name of the timer.
 *
 * @return The timer.
 */
enum spamm_timer_type_t
spamm_timer_get_timer_type (const char *name)
{
  enum spamm_timer_type_t timer_type = -1;

  if (strcasecmp(name, "walltime") == 0)
  {
    timer_type = walltime;
  }

  else if (strcasecmp(name, "papi_total_cycles") == 0)
  {
    timer_type = papi_total_cycles;
  }

  else if (strcasecmp(name, "papi_flop") == 0)
  {
    timer_type = papi_flop;
  }

  else if (strcasecmp(name, "papi_l1_dcm") == 0)
  {
    timer_type = papi_l1_dcm;
  }

  else if (strcasecmp(name, "papi_l1_icm") == 0)
  {
    timer_type = papi_l1_icm;
  }

  else if (strcasecmp(name, "papi_l2_dcm") == 0)
  {
    timer_type = papi_l2_dcm;
  }

  else if (strcasecmp(name, "papi_l2_icm") == 0)
  {
    timer_type = papi_l2_icm;
  }

  else
  {
    timer_type = -1;
  }

  return timer_type;
}

/** Get native events and print a list.
 */
void
spamm_timer_get_native_events ()
{
#ifdef HAVE_PAPI
  int i;
  int event_number = 0;
  int retval;
  PAPI_event_info_t info;

  /* Initialize PAPI library. */
  spamm_timer_PAPI_init();

  printf("Native events:\n");
  i = PAPI_NATIVE_MASK;
  printf("Nr.   Name                           Code        Description\n");
  do
  {
    if ((retval = PAPI_get_event_info(i, &info)) == PAPI_OK)
    {
      printf("%3u   %-30s 0x%-10x%s\n", event_number++, info.symbol, info.event_code, info.long_descr);
    }
  }
  while ((retval = PAPI_enum_event(&i, PAPI_ENUM_ALL)) == PAPI_OK);

  printf("Preset events:\n");
  i = PAPI_PRESET_MASK;
  printf("Nr.   Name                           Code        Description\n");
  do
  {
    if ((retval = PAPI_get_event_info(i, &info)) == PAPI_OK)
    {
      printf("%3u   %-30s 0x%-10x%s\n", event_number++, info.symbol, info.event_code, info.long_descr);
    }
  }
  while ((retval = PAPI_enum_event(&i, PAPI_ENUM_ALL)) == PAPI_OK);
#else
  printf("PAPI is not available\n");
  exit(1);
#endif
}
