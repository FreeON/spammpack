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

  /** The length of the values array. */
  short length;

  /** An array that holds the names of the PAPI events. */
  char **event_names;

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
  spamm_timer_PAPI_init();
#endif

  switch (timer->type)
  {
    case walltime:
    case empty:
      break;

#ifdef HAVE_PAPI
    case papi_total_instructions:
    case papi_total_cycles:
    case papi_flop:
    case papi_vec_sp:
    case papi_l1_dcm:
    case papi_l1_icm:
    case papi_l2_dcm:
    case papi_l2_icm:
      if ((papi_result = PAPI_create_eventset(&timer->eventset)) != PAPI_OK)
      {
        spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_create_eventset()");
      }

      switch(timer->type)
      {
        case walltime:
          /* We already took care of this option, just making the compiler
           * happy.
           */
          break;

        case papi_total_instructions:
          if ((papi_result = PAPI_add_event(timer->eventset, PAPI_TOT_INS)) != PAPI_OK)
          {
            spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_add_event()");
          }

          timer->length = 1;
          timer->event_values = calloc(timer->length, sizeof(long_long));
          timer->event_names = calloc(timer->length, sizeof(char*));
          timer->event_names[0] = strndup("PAPI_TOT_INS", strlen("PAPI_TOT_INS"));
          break;

        case papi_total_cycles:
          if ((papi_result = PAPI_add_event(timer->eventset, PAPI_TOT_CYC)) != PAPI_OK)
          {
            spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_add_event()");
          }

          timer->length = 1;
          timer->event_values = calloc(timer->length, sizeof(long_long));
          timer->event_names = calloc(timer->length, sizeof(char*));
          timer->event_names[0] = strndup("PAPI_TOT_CYC", strlen("PAPI_TOT_CYC"));
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

          timer->length = 2;
          timer->event_values = calloc(timer->length, sizeof(long_long));
          timer->event_names = calloc(timer->length, sizeof(char*));
          timer->event_names[0] = strndup("PAPI_FP_OPS", strlen("PAPI_FP_OPS"));
          timer->event_names[0] = strndup("PAPI_TOT_CYC", strlen("PAPI_TOT_CYC"));
          break;

        case papi_vec_sp:
          if ((papi_result = PAPI_add_event(timer->eventset, PAPI_VEC_SP)) != PAPI_OK)
          {
            spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_add_event()");
          }

          timer->length = 1;
          timer->event_values = calloc(timer->length, sizeof(long_long));
          timer->event_names = calloc(timer->length, sizeof(char*));
          timer->event_names[0] = strndup("PAPI_VEC_SP", strlen("PAPI_VEC_SP"));
          break;

        case papi_l1_dcm:
          if ((papi_result = PAPI_add_event(timer->eventset, PAPI_L1_DCM)) != PAPI_OK)
          {
            spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_add_event()");
          }

          timer->length = 1;
          timer->event_values = calloc(timer->length, sizeof(long_long));
          timer->event_names = calloc(timer->length, sizeof(char*));
          timer->event_names[0] = strndup("PAPI_L1_DCM", strlen("PAPI_L1_DCM"));
          break;

        case papi_l1_icm:
          if ((papi_result = PAPI_add_event(timer->eventset, PAPI_L1_ICM)) != PAPI_OK)
          {
            spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_add_event()");
          }

          timer->length = 1;
          timer->event_values = calloc(timer->length, sizeof(long_long));
          timer->event_names = calloc(timer->length, sizeof(char*));
          timer->event_names[0] = strndup("PAPI_L1_ICM", strlen("PAPI_L1_ICM"));
          break;

        case papi_l2_dcm:
          if ((papi_result = PAPI_add_event(timer->eventset, PAPI_L2_DCM)) != PAPI_OK)
          {
            spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_add_event()");
          }

          timer->length = 1;
          timer->event_values = calloc(timer->length, sizeof(long_long));
          timer->event_names = calloc(timer->length, sizeof(char*));
          timer->event_names[0] = strndup("PAPI_L2_DCM", strlen("PAPI_L2_DCM"));
          break;

        case papi_l2_icm:
          if ((papi_result = PAPI_add_event(timer->eventset, PAPI_L2_ICM)) != PAPI_OK)
          {
            spamm_timer_handle_PAPI_error(papi_result, "new timer, PAPI_add_event()");
          }

          timer->length = 1;
          timer->event_values = calloc(timer->length, sizeof(long_long));
          timer->event_names = calloc(timer->length, sizeof(char*));
          timer->event_names[0] = strndup("PAPI_L2_ICM", strlen("PAPI_L2_ICM"));
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
  int i;
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
    case papi_l1_dcm:
    case papi_l1_icm:
    case papi_l2_dcm:
    case papi_l2_icm:
      if ((papi_result = PAPI_cleanup_eventset((*timer)->eventset)) != PAPI_OK)
      {
        spamm_timer_handle_PAPI_error(papi_result, "delete timer, PAPI_cleanup_eventset()");
      }

      if ((papi_result = PAPI_destroy_eventset(&(*timer)->eventset)) != PAPI_OK)
      {
        spamm_timer_handle_PAPI_error(papi_result, "delete timer, PAPI_destroy_eventset()");
      }

      for (i = 0; i < (*timer)->length; i++)
      {
        free((*timer)->event_names[i]);
      }
      free((*timer)->event_names);
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
    case papi_l1_dcm:
    case papi_l1_icm:
    case papi_l2_dcm:
    case papi_l2_icm:
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
    case papi_l1_dcm:
    case papi_l1_icm:
    case papi_l2_dcm:
    case papi_l2_icm:
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
  assert(length != NULL);
  assert(values != NULL);

  if (timer->timer_running != 0)
  {
    printf("[timer get] this timer is running right now\n");
    exit(1);
  }

  switch (timer->type)
  {
    case walltime:
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
      break;

#ifdef HAVE_PAPI
    case papi_total_instructions:
    case papi_total_cycles:
    case papi_vec_sp:
    case papi_l1_dcm:
    case papi_l1_icm:
    case papi_l2_dcm:
    case papi_l2_icm:
      *length = 1;
      *values = calloc(*length, sizeof(unsigned long long));
      (*values)[0] = timer->event_values[0];
      break;

    case papi_flop:
      *length = 2;
      *values = calloc(*length, sizeof(unsigned long long));
      (*values)[0] = timer->event_values[0];
      (*values)[1] = timer->event_values[1];
      break;
#endif

    default:
      printf("[timer get] unknown timer type\n");
      exit(1);
      break;
  }
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
  int i;
  char *result = calloc(2000, sizeof(char));

#ifdef HAVE_PAPI
  sprintf(result, "[ ");
  for (i = 0; i < timer->length; i++)
  {
    sprintf(result, "%s %lli (%s)", result, timer->event_values[i], timer->event_names[i]);
  }
  sprintf(result, " ]");
#else
  sprintf(result, "walltime %u", 0);
#endif

  return result;
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

    case papi_l1_dcm:
      sprintf(string, "PAPI Level 1 data cache misses");
      break;

    case papi_l1_icm:
      sprintf(string, "PAPI Level 1 instruction cache misses");
      break;

    case papi_l2_dcm:
      sprintf(string, "PAPI Level 2 data cache misses");
      break;

    case papi_l2_icm:
      sprintf(string, "PAPI Level 2 instruction cache misses");
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
