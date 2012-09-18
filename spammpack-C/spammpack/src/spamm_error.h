/** @file */

#ifndef __SPAMM_ERROR_H
#define __SPAMM_ERROR_H

/** A short-hand for a fatal error message. */
#define SPAMM_FATAL(...) spamm_error_fatal(__FILE__, __LINE__, __VA_ARGS__)

/** A short-hand for a warning message. */
#define SPAMM_WARN(...) spamm_error_warning(__FILE__, __LINE__, __VA_ARGS__)

void
spamm_error_fatal (const char *const filename, const int line, ...);

void
spamm_error_warning (const char *const filename, const int line, ...);

#endif
