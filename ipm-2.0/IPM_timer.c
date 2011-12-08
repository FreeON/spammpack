const char* IPM_timer_rcsid = "$Header: /u/spert/src/IPM/RCS/IPM_timer.c,v 1.37 1996/08/08 00:43:34 krste Exp $";
/*
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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef EXIT_SUCCESS
# define EXIT_SUCCESS (0)
#endif
#ifndef EXIT_FAILURE
# define EXIT_FAILURE (0)
#endif

/******************************************************************/
/* Configure must set up these symbols depending on system type.  */
/******************************************************************/

/* Must define at least one module. */
#if defined(HAVE_GETRUSAGE)
#elif defined(HAVE_GETTIMEOFDAY)
#elif defined(HAVE_SYSSGI)
#elif defined(HAVE_GETHRTIME)
#elif defined(HAVE_GETTIMER)
#elif defined(HAVE_T0VCOUNT)
#else
error
#endif

/* Must have defined the default module. */

#ifndef IPM_init_DEFAULT
error
#endif

/* Include IPM header files after configured symbols. */

#include "IPM_timer.h"
#include "IPM_timer_internal.h"

#ifdef HAVE_GETRUSAGE
#include "IPM_timer_getrusage.h"
#endif

#ifdef HAVE_GETTIMEOFDAY
#include "IPM_timer_gettimeofday.h"
#endif

#ifdef HAVE_T0VCOUNT
#include "IPM_timer_t0vcount.h"
#endif

#ifdef HAVE_SYSSGI
#include "IPM_timer_syssgi.h"
#endif

#ifdef HAVE_GETHRTIME
#include "IPM_timer_gethrtime.h"
#endif

#ifdef HAVE_GETTIMER
#include "IPM_timer_gettimer.h"
#endif

/*

    Common functions to handle resolution and range reporting.
    
    By default, modules just initialize these global variables.
    Modules can install their own handlers if more support required.

*/    

double IPM_timer_resolution_val;
double IPM_timer_range_val;

static double
IPM_timer_resolution_common(void)
{ return IPM_timer_resolution_val; }

static double
IPM_timer_range_common(void)
{ return IPM_timer_range_val; }

/*

    Functions that handle timer format string.  Format string is should not
    overflow with largest time, and also round to available precision.
    Always prints times in format iiiiii.ffffff where i is integer seconds,
    and f is fractional seconds.

*/

/* Space for format string. */
#define FORMAT_LEN (32)
static char IPM_timer_format_str[FORMAT_LEN];

static void
IPM_timer_format_set(double resolution, double range)
{
    /*
    Simple minded approach to finding ranges.  Don't want to have to link
    with math library to get log10 functions.
    */

    int int_digits = 1;         /* Always at least one integer digit. */
    int frac_digits = 0;
    double magnitude = 10.0;

    assert(range > 0 && range < 1e15);

    while (magnitude <= range)
    {
        magnitude *= 10.0;
        int_digits++;
    }

    assert(resolution >= 1e-12); /* Nothing better than picosecond
                                    resolution. */

    while (resolution <= 0.75)
    {
        resolution *= 10.0;
        frac_digits++;
    }

    (void) sprintf(IPM_timer_format_str,
                            "%%%d.%df",
                            int_digits + frac_digits + 1,
                            frac_digits);
}

static char*
IPM_timer_format_common(void)
{
    return IPM_timer_format_str;
}

/*

    Common function to handle timer reports.

*/

char* IPM_timer_module = NULL;
char* IPM_timer_type = NULL;

static void
IPM_timer_report_common(const char* tag, const IPM_timer* t, char* buf)
{
    /* Uses strlen() hack to get around sprintf not always returning
       character count. */

    char* bufp = buf;

    (void) sprintf(bufp, "IPM_timer: ");
    bufp = buf + strlen(buf);

    if (tag!=NULL)
    {
        (void) sprintf(bufp, "tag=%s ", tag);
        bufp = buf + strlen(buf);
    }

    (void) sprintf(bufp, "module=%s type=%s time=",
                   IPM_timer_module, IPM_timer_type);
    bufp = buf + strlen(buf);

    (void) sprintf(bufp, IPM_timer_format(), IPM_timer_read(t));
    bufp = buf + strlen(buf);

    (void) sprintf(bufp, "\n");
}

/* Standard error reporting functions. */

void
IPM_error(const char *message)
{
    fprintf(stderr, "IPM %s: ERROR: ", IPM_version);

    if (IPM_timer_module!=NULL)
        fprintf(stderr, "module=%s: ", IPM_timer_module);

    if (IPM_timer_type!=NULL)
        fprintf(stderr, "type=%s: ", IPM_timer_type);

    fprintf(stderr, "%s\n", message);

    exit(EXIT_FAILURE);
}

void
IPM_system_error(const char*message)
{
    perror("IPM: SYSTEM ERROR: ");
    IPM_error(message);
}

/*

    Initialization function parses environment, sets up selected timing
    module in dispatch table, and calls module's initialization code.

*/

static void
IPM_init(void)
{
    static int already_initialized = 0;
    char *mp;

    /* If we execute this code twice, the timer module must not have
       initialized all the function pointers correctly. */
    assert(already_initialized==0);
    already_initialized = 1;

    /* Install default handlers. */

    IPM_timer_resolution_p = IPM_timer_resolution_common;
    IPM_timer_range_p = IPM_timer_range_common;
    IPM_timer_format_p = IPM_timer_format_common;
    IPM_timer_report_p = IPM_timer_report_common;

    mp = getenv("IPM_timer_module");

    IPM_timer_module = mp;

    if (!mp)
        IPM_init_DEFAULT();

#ifdef HAVE_GETRUSAGE
    else if (!strcmp(mp, "getrusage"))
        IPM_init_getrusage();
#endif

#ifdef HAVE_GETTIMEOFDAY
    else if (!strcmp(mp, "gettimeofday"))
        IPM_init_gettimeofday();
#endif

#ifdef HAVE_T0VCOUNT
    else if (!strcmp(mp, "t0vcount"))
        IPM_init_t0vcount();
#endif

#ifdef HAVE_SYSSGI
    else if (!strcmp(mp, "syssgi"))
        IPM_init_syssgi();
#endif

#ifdef HAVE_GETHRTIME
    else if (!strcmp(mp, "gethrtime"))
        IPM_init_gethrtime();
#endif

#ifdef HAVE_GETTIMER
    else if (!strcmp(mp, "gettimer"))
        IPM_init_gettimer();
#endif

    else
    {
        IPM_error("No such IPM_timer_module.");
    }

    /* Set up format string. */
    IPM_timer_format_set(IPM_timer_resolution(), IPM_timer_range());

    /* Check these were set by init routine. */
    assert(IPM_timer_module!=NULL);
    assert(IPM_timer_type!=NULL);
}

static double
IPM_timer_resolution_init(void)
{
    IPM_init();                    /* Set things up, */
    return IPM_timer_resolution(); /* then try again. */
}

static double
IPM_timer_range_init(void)
{
    IPM_init();               /* Set things up, */
    return IPM_timer_range(); /* then try again. */
}


static double
IPM_timer_read_init(const IPM_timer* t)
{
    IPM_error("IPM_timer_read() called before measurement.");
    return 0.0;             /* Keep compiler happy. */
}

static char*
IPM_timer_format_init(void)
{
    IPM_init();                 /* First call, so set up everything. */
    return IPM_timer_format();  /* Now return value. */
}

static void
IPM_timer_report_init(const char* tag, const IPM_timer* t, char* buf)
{
    IPM_error("IPM_timer_report() called before measurement.");
}

static void
IPM_timer_clear_init(IPM_timer* t)
{
    IPM_init();             /* Set things up, */
    IPM_timer_clear(t);     /* then try again. */
}

static void
IPM_timer_start_init(IPM_timer* t)
{
    IPM_init();             /* Set things up, */
    IPM_timer_start(t);     /* then try again. */
}

static void
IPM_timer_stop_init(IPM_timer* t)
{
    IPM_error("Called IPM_stop() before IPM_start().");
}

static double
IPM_timer_get_time_init(void)
{
    IPM_init();                  /* Set things up, */
    return IPM_timer_get_time(); /* then try again. */
}

/*
    Initialize the dispatch table with the *_init functions.
*/

IPM_timer_resolution_p_t IPM_timer_resolution_p = IPM_timer_resolution_init;
IPM_timer_range_p_t IPM_timer_range_p = IPM_timer_range_init;
IPM_timer_clear_p_t IPM_timer_clear_p = IPM_timer_clear_init;
IPM_timer_start_p_t IPM_timer_start_p = IPM_timer_start_init;
IPM_timer_stop_p_t IPM_timer_stop_p = IPM_timer_stop_init;
IPM_timer_read_p_t IPM_timer_read_p = IPM_timer_read_init;
IPM_timer_format_p_t IPM_timer_format_p = IPM_timer_format_init;
IPM_timer_report_p_t IPM_timer_report_p = IPM_timer_report_init;
IPM_timer_get_time_p_t IPM_timer_get_time_p = IPM_timer_get_time_init;
