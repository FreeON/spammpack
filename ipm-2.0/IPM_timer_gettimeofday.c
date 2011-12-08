const char* IPM_timer_gettimeofday_rcsid = "$Header: /u/spert/src/IPM/RCS/IPM_timer_gettimeofday.c,v 1.18 1996/08/08 00:43:40 krste Exp $";

/*
    Code that uses gettimeofday() to implement an IPM timer.
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

#include <sys/time.h>
#include <sys/resource.h>

#include "IPM_timer.h"
#include "IPM_timer_internal.h"
#include "IPM_timer_gettimeofday.h"

/*****************************************************************/
/* Configure can set up these symbols depending on system type.  */
/*****************************************************************/

/* Most systems have microsecond resolution. Configure changes value if
   not. */
#ifndef IPM_GETTIMEOFDAY_RESOLUTION
#  define IPM_GETTIMEOFDAY_RESOLUTION (1e-6)
#endif

/* 2^31-1 is maximum number of +ve seconds in a 32-bit timeval. */
#ifndef IPM_GETTIMEOFDAY_RANGE
# define IPM_GETTIMEOFDAY_RANGE (2147483647.0)
#endif

typedef struct
{
    struct timeval start;       /* Start values. */
    struct timeval stop;        /* Stop values. */
    double t_us;                /* Accumulated microseconds, stored in */
                                /* double. */
} IPM_timer_gettimeofday;

static void
IPM_error_gettimeofday(const char* message)
{
    perror(NULL);
    IPM_error(message);
}

static void
IPM_timer_clear_gettimeofday(IPM_timer *ipmt)
{
    IPM_timer_gettimeofday *gtodt = (IPM_timer_gettimeofday *) ipmt;

    gtodt->t_us = 0.0;              /* Clear accumulated time. */
}

static void
IPM_timer_start_gettimeofday(IPM_timer *ipmt)
{
    IPM_timer_gettimeofday *gtodt = (IPM_timer_gettimeofday *) ipmt;

    if (gettimeofday(&(gtodt->start), NULL))
        IPM_error_gettimeofday("IPM_timer_start()");
}

/* Do timeval calculations in microseconds as this avoids the accuracy loss */
/* that would occur when converting 1e-6 to binary if working in seconds. */
#define TV2US(timeval) (((timeval).tv_sec * 1e6) + (timeval).tv_usec)

static void
IPM_timer_stop_gettimeofday(IPM_timer *ipmt)
{
    IPM_timer_gettimeofday *gtodt = (IPM_timer_gettimeofday *) ipmt;
  
    if (gettimeofday(&(gtodt->stop), NULL))
        IPM_error_gettimeofday("IPM_timer_stop()");

    gtodt->t_us += TV2US(gtodt->stop) - TV2US(gtodt->start);
}

static double
IPM_timer_read_gettimeofday(const IPM_timer *ipmt)
{
    const IPM_timer_gettimeofday *gtodt
        = (const IPM_timer_gettimeofday *) ipmt;

    return gtodt->t_us * 1e-6;
}

static double
IPM_timer_get_time_gettimeofday(void)
{
    struct timeval t;

    if (gettimeofday(&t, NULL))
        IPM_error_gettimeofday("IPM_timer_get_time()");

    return TV2US(t) * 1e-6;
}

void
IPM_init_gettimeofday(void)
{
    /* Make sure opaque struct large enough. */
    assert(sizeof(IPM_timer_gettimeofday) <= IPM_TIMER_STRUCT_MAX_SIZE);

    IPM_timer_module = "gettimeofday";
    IPM_timer_type = "real";

    /* Initialize resolution and range variables. */
    IPM_timer_resolution_val = IPM_GETTIMEOFDAY_RESOLUTION;
    IPM_timer_range_val = IPM_GETTIMEOFDAY_RANGE;

    /* Install function pointers. */
    IPM_timer_clear_p = IPM_timer_clear_gettimeofday;
    IPM_timer_start_p = IPM_timer_start_gettimeofday;
    IPM_timer_stop_p = IPM_timer_stop_gettimeofday;
    IPM_timer_read_p = IPM_timer_read_gettimeofday;
    IPM_timer_get_time_p = IPM_timer_get_time_gettimeofday;
}
