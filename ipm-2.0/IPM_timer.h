/*
    $Header: /u/spert/src/IPM/RCS/IPM_timer.h,v 1.29 1996/08/07 21:13:19 krste Exp $

    IPM timer package provides portable, high resolution interval timing.
*/

/*
    "Copyright (c) 1995 The Regents of the University of California.  
    All rights reserved."

    Permission to use, copy, modify, and distribute this software and
    its documentation for any purpose, without fee, and without
    written agreement is hereby granted, provided that the above
    copyright notice and the following two paragraphs appear in all
    copies of this software.

    IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY
    PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL
    DAMAGES ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
    DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN
    ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY
    WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
    OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE
    SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE
    UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE,
    SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS."

    Author: Krste Asanovic <krste@cs.berkeley.edu>
*/

#ifndef IPM_timer_h_
#define IPM_timer_h_

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque structure defined here. */
#define IPM_TIMER_STRUCT_MAX_SIZE (1024)

/*
  HACK: Makes assumption that double has worst alignment.  Don't want to
  use long double as not all compilers support it.
*/
typedef struct {
    union {
        double worst_alignment;
        char storage[IPM_TIMER_STRUCT_MAX_SIZE];
    } timer_union;
} IPM_timer;

/* Declare function dispatch table. */
typedef double (*IPM_timer_resolution_p_t)(void);
typedef double (*IPM_timer_range_p_t)(void);
typedef void (*IPM_timer_clear_p_t)(IPM_timer *);
typedef void (*IPM_timer_start_p_t)(IPM_timer *);
typedef void (*IPM_timer_stop_p_t)(IPM_timer *);
typedef double (*IPM_timer_read_p_t)(const IPM_timer *);
typedef char* (*IPM_timer_format_p_t)(void);
typedef void (*IPM_timer_report_p_t)(const char*, const IPM_timer *, char*);
typedef double (*IPM_timer_get_time_p_t)(void);

extern IPM_timer_resolution_p_t IPM_timer_resolution_p;
extern IPM_timer_range_p_t IPM_timer_range_p;
extern IPM_timer_clear_p_t IPM_timer_clear_p;
extern IPM_timer_start_p_t IPM_timer_start_p;
extern IPM_timer_stop_p_t IPM_timer_stop_p;
extern IPM_timer_read_p_t IPM_timer_read_p;
extern IPM_timer_format_p_t IPM_timer_format_p;
extern IPM_timer_report_p_t IPM_timer_report_p;
extern IPM_timer_get_time_p_t IPM_timer_get_time_p;

/* Declare macros to call functions through dispatch table. */

/* Resolution of timer clock given in seconds per clock tick. */
#define IPM_timer_resolution() ((*IPM_timer_resolution_p)())
       
/* Maximum interval in seconds that can be measured before timer overflows. */
#define IPM_timer_range() ((*IPM_timer_range_p)())

/* Reset accumulated time in interval timer. */
#define IPM_timer_clear(t) ((*IPM_timer_clear_p)(t))

/* Start timer. */
#define IPM_timer_start(t) ((*IPM_timer_start_p)(t))

/* Stop timer, add on time to accumulated time. */
#define IPM_timer_stop(t) ((*IPM_timer_stop_p)(t))

/* Return accumulated interval time in seconds. */
#define IPM_timer_read(t) ((*IPM_timer_read_p)(t))

/* Return format string appropriate for printing timer values. */
#define IPM_timer_format() ((*IPM_timer_format_p)())

/*
    IPM_timer_report(const char *tag, const IPM_timer *t, char *buf);

    Returns string given timer value in buf. The string buf must be at
    least IPM_TIMER_REPORT_LEN + strlen(tag) bytes long. The tag is a name
    that is used to identify the timer report.
*/

#define IPM_TIMER_REPORT_LEN (200)

#define IPM_timer_report(tag, t, buf) ((*IPM_timer_report_p)(tag, t, buf))

/* Return current timer clock time in seconds. */
#define IPM_timer_get_time() (*IPM_timer_get_time_p)()

#ifdef __cplusplus
}
#endif

/* End of #ifndef IPM_timer_h_ */
#endif
