#include "spamm.h"
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** The global loglevel. */
static enum spamm_log_severity_t spamm_loglevel = fatal;

/** Set the global loglevel.
 *
 * @param loglevel The loglevel to set.
 */
void
spamm_set_loglevel (const enum spamm_log_severity_t loglevel)
{
  if (loglevel >= 0 && loglevel <= debug)
  {
    spamm_loglevel = loglevel;
  }

  else
  {
    LOG_FATAL("illegal value for loglevel (%i)\n", loglevel);
    exit(1);
  }
}

/** Logging function.
 *
 * Use spamm_log() to print out messages. The messages are prepended with
 * filename and linenumber. Typical use:
 *
 * <code>spamm_log("opening new file: %s", __FILE__, __LINE__, filename);</code>
 *
 * @param severity The severity of this message.
 * @param format Format string. See printf() for a detailed description of its
 *               syntax.
 * @param filename The filename of the source of the caller.
 * @param linenumber The linenumber of the caller.
 */
void
spamm_log (const enum spamm_log_severity_t severity, const char *format,
    const char *filename, const unsigned int linenumber, ...)
{
  char *new_format;
  char *severity_string[] = { "FATAL", "INFO", "DEBUG" };
  int size;
  va_list ap;

  assert(format != NULL);
  assert(filename != NULL);

  if (severity <= spamm_loglevel)
  {
    va_start(ap, linenumber);

    /* Fix up the format, i.e. insert file and line number. We conservatively
     * estimate that the file name and the line number will not exceed 2000
     * characters. */
    size = sizeof(char)*(strlen(format)+2000);
    new_format = (char*) malloc(size);
    snprintf(new_format, size, "[%s:%u - %s] %s", filename, linenumber, severity_string[severity], format);

    vprintf(new_format, ap);
    fflush(stdout);

    free(new_format);
    va_end(ap);
  }
}
