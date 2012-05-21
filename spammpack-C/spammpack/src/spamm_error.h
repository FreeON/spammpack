/** @file */

#ifndef __SPAMM_ERROR_H
#define __SPAMM_ERROR_H

#define SPAMM_FATAL(format, ...) spamm_error_fatal(__FILE__, __LINE__, format, __VA_ARGS__)

void
spamm_error_fatal (const char *const filename, const int line, const char *const format, ...);

#endif
