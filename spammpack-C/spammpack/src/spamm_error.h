/** @file */

#ifndef __SPAMM_ERROR_H
#define __SPAMM_ERROR_H

/** A short-hand for a fatal error message. */
#define SPAMM_FATAL(...) spamm_error_fatal(__FILE__, __LINE__, __VA_ARGS__)

/** A short-hand for a warning message. */
#define SPAMM_WARN(...) spamm_error_warning(__FILE__, __LINE__, __VA_ARGS__)

/** A short-hand for an info message. */
#define SPAMM_INFO(...) spamm_error_info(__FILE__, __LINE__, __VA_ARGS__)

#ifdef __cplusplus
#define __BEGIN_DECLARATIONS extern "C" {
#define __END_DECLARATIONS }
#else
#define __BEGIN_DECLARATIONS
#define __END_DECLARATIONS
#endif

__BEGIN_DECLARATIONS

void
spamm_error_fatal (const char *const filename, const int line, ...);

void
spamm_error_warning (const char *const filename, const int line, ...);

void
spamm_error_info (const char *const filename, const int line, ...);

__END_DECLARATIONS

#endif
