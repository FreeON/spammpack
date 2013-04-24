#include "spamm_error.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <execinfo.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
static short print_lock_is_initialized = 0;
static omp_lock_t spamm_print_lock;
#endif

/** Print a backtrace. */
void
spamm_error_print_backtrace ()
{
  void *array[50];
  size_t size;
  char **strings;
  size_t i;

  size = backtrace(array, 50);
  strings = backtrace_symbols(array, size);

  printf("Obtained %zd stack frames.\n", size-2);

  /* Omit the first 2 frames since we know those are the error handlers. */
  for(i = 2; i < size; i++)
  {
    printf("%s\n", strings[i]);
  }

  free(strings);
}

/** Print out a warning message.
 *
 * @param filename The filename this error message was generated in.
 * @param line The line number this error message was generated on.
 * @param format The format of the error message.
 */
void
spamm_error_warning (const char *const filename, const int line, ...)
{
  va_list va;
  int format_length;

  char *old_format;
  char *new_format;

#ifdef _OPENMP
  /* Initialize lock. */
  if(!print_lock_is_initialized)
  {
    omp_init_lock(&spamm_print_lock);
    print_lock_is_initialized = 1;
  }
#endif

  /* Initialize variadic argument list. */
  va_start(va, line);

  /* Get format string. */
  old_format = va_arg(va, char*);

  format_length = 6+strlen(filename)+100+strlen(old_format);
  new_format = calloc(format_length, sizeof(char));
#ifdef _OPENMP
  snprintf(new_format, format_length-1, "[%s:%i thread %i WARNING] ", filename, line, omp_get_thread_num());
#else
  snprintf(new_format, format_length-1, "[%s:%i WARNING] ", filename, line);
#endif
  strncat(new_format, old_format, format_length);

  /* Print error. */
#ifdef _OPENMP
  omp_set_lock(&spamm_print_lock);
#endif
  vprintf(new_format, va);
#ifdef _OPENMP
  omp_unset_lock(&spamm_print_lock);
#endif

  /* Cleanup. */
  va_end(va);
  free(new_format);
}

/** Print out a fatal error message.
 *
 * @param filename The filename this error message was generated in.
 * @param line The line number this error message was generated on.
 * @param format The format of the error message.
 */
void
spamm_error_fatal (const char *const filename, const int line, ...)
{
  va_list va;
  int format_length;

  char *old_format;
  char *new_format;

#ifdef _OPENMP
  /* Initialize lock. */
  if(!print_lock_is_initialized)
  {
    omp_init_lock(&spamm_print_lock);
    print_lock_is_initialized = 1;
  }
#endif

  /* Initialize variadic argument list. */
  va_start(va, line);

  /* Get format string. */
  old_format = va_arg(va, char*);

  format_length = 6+strlen(filename)+100+strlen(old_format);
  new_format = calloc(format_length, sizeof(char));
#ifdef _OPENMP
  snprintf(new_format, format_length-1, "[%s:%i thread %i FATAL] ", filename, line, omp_get_thread_num());
#else
  snprintf(new_format, format_length-1, "[%s:%i FATAL] ", filename, line);
#endif
  strncat(new_format, old_format, format_length);

  /* Print error. */
#ifdef _OPENMP
  omp_set_lock(&spamm_print_lock);
#endif
  vprintf(new_format, va);
#ifdef _OPENMP
  omp_unset_lock(&spamm_print_lock);
#endif

  /* Cleanup. */
  va_end(va);
  free(new_format);

  /* Print backtrace. */
  printf("\n");
  spamm_error_print_backtrace();

  /* Exit with signal so debuggers can produce a backtrace. */
  abort();
}
