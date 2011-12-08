const char* IPM_timer_gethrtime_rcsid = "$Header: /u/spert/src/IPM/RCS/IPM_timer_gethrtime.c,v 1.25 1996/08/08 00:43:36 krste Exp $";

/*
    Code that uses the gethrtime()/gethrvtime() calls on SunOS5 systems to
    implement an IPM timer.
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

/********************************************************/
/* Configure can set up these system dependent values.  */
/********************************************************/

#ifndef IPM_GETHRTIME_RANGE
   /* 64-bit counter with 1 ns resolution. */
#  define IPM_GETHRTIME_RANGE (18446744073.6)
#endif

#ifndef IPM_GETHRTIME_RESOLUTION
   /* Apparently, this is the resolution on most Sun systems. */
#  define IPM_GETHRTIME_RESOLUTION (500e-9)
#endif

#ifndef IPM_GETHRVTIME_RESOLUTION
   /* 100 Hz quanta on SunOS 5 SPARC systems. */
#  define IPM_GETHRVTIME_RESOLUTION (0.010)
#endif

#include "IPM_timer.h"
#include "IPM_timer_internal.h"
#include "IPM_timer_gethrtime.h"

typedef struct
{
    hrtime_t start;
    hrtime_t total;
} IPM_timer_gethrtime;

static void
IPM_timer_clear_gethrtime(IPM_timer *ipmt)
{
    IPM_timer_gethrtime *ghrt = (IPM_timer_gethrtime *) ipmt;
    
    ghrt->total = 0;
}

static void
IPM_timer_start_gethrtime(IPM_timer *ipmt)
{
    IPM_timer_gethrtime *ghrt = (IPM_timer_gethrtime *) ipmt;

    ghrt->start = gethrtime();
}

static void
IPM_timer_stop_gethrtime(IPM_timer *ipmt)
{
    IPM_timer_gethrtime *ghrt = (IPM_timer_gethrtime *) ipmt;

    const hrtime_t stop = gethrtime();
    
    ghrt->total += (stop - ghrt->start);
}

static void
IPM_timer_start_gethrvtime(IPM_timer *ipmt)
{
    IPM_timer_gethrtime *ghrt = (IPM_timer_gethrtime *) ipmt;

    ghrt->start = gethrvtime();
}

static void
IPM_timer_stop_gethrvtime(IPM_timer *ipmt)
{
    IPM_timer_gethrtime *ghrt = (IPM_timer_gethrtime *) ipmt;

    const hrtime_t stop = gethrvtime();
    
    ghrt->total += (stop - ghrt->start);
}

static double
IPM_timer_read_gethrtime(const IPM_timer *ipmt)
{
    const IPM_timer_gethrtime *ghrt = (const IPM_timer_gethrtime *) ipmt;

    return ghrt->total * 1e-9;
}

static double
IPM_timer_get_time_gethrtime(void)
{
    hrtime_t t;

    t = gethrtime();

    return t * 1e-9;
}

static double
IPM_timer_get_time_gethrvtime(void)
{
    hrtime_t t;

    t = gethrvtime();

    return t * 1e-9;
}

void
IPM_init_gethrtime(void)
{
    char *timer_str;

    /* Make sure opaque struct large enough. */
    assert(sizeof(IPM_timer_gethrtime) <= IPM_TIMER_STRUCT_MAX_SIZE);

    IPM_timer_module = "gethrtime";
    IPM_timer_range_val = IPM_GETHRTIME_RANGE;

    /* Parse environment variables. */
    timer_str = getenv("IPM_timer_type");
    IPM_timer_type = timer_str;

    if (!timer_str || !strcmp(timer_str, "real"))
    {
        IPM_timer_resolution_val = IPM_GETHRTIME_RESOLUTION;

        IPM_timer_type = "real";
        IPM_timer_start_p = IPM_timer_start_gethrtime;
        IPM_timer_stop_p = IPM_timer_stop_gethrtime;
        IPM_timer_get_time_p = IPM_timer_get_time_gethrtime;
    }
    else if (!strcmp(timer_str, "process"))
    {
        IPM_timer_resolution_val = IPM_GETHRVTIME_RESOLUTION;

        IPM_timer_type = "process";
        IPM_timer_start_p = IPM_timer_start_gethrvtime;
        IPM_timer_stop_p = IPM_timer_stop_gethrvtime;
        IPM_timer_get_time_p = IPM_timer_get_time_gethrvtime;
    }
    else
        IPM_error("IPM_init() : No such IPM_timer_type.");

    /* Install remaining function pointers. */
    IPM_timer_clear_p = IPM_timer_clear_gethrtime;
    IPM_timer_read_p = IPM_timer_read_gethrtime;
}
