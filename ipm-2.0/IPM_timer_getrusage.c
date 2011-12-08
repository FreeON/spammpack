const char* IPM_timer_getrusage_rcsid = "$Header: /u/spert/src/IPM/RCS/IPM_timer_getrusage.c,v 1.23 1996/08/08 00:43:38 krste Exp $";

/*
    Code that uses BSD getrusage() to implement an IPM timer.
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "IPM_timer.h"
#include "IPM_timer_internal.h"
#include "IPM_timer_getrusage.h"

/******************************************************************/
/* Configure should set up these symbols depending on system type.  */
/******************************************************************/

/* Most getrusages have 100Hz clock. */
#ifndef IPM_GETRUSAGE_RESOLUTION
#  define IPM_GETRUSAGE_RESOLUTION (0.010)
#endif

/* 2^31-1 is maximum number of +ve seconds in 32-bit timeval. */
#ifndef IPM_GETRUSAGE_RANGE
#  define IPM_GETRUSAGE_RANGE (2147483647.0)
#endif

typedef struct
{
    struct rusage start;        /* Start values. */
    struct rusage stop;         /* Stop values. */
    double t_us;                /* Accumulated microseconds, stored in */
                                /* double. */
} IPM_timer_getrusage;

static void
IPM_error_getrusage(const char* message)
{
    perror(NULL);
    IPM_error(message);
}

static void
IPM_timer_clear_getrusage(IPM_timer* ipmt)
{
    IPM_timer_getrusage *grut = (IPM_timer_getrusage *) ipmt;

    grut->t_us = 0.0;              /* Clear accumulated time. */
}

static void
IPM_timer_start_getrusage(IPM_timer* ipmt)
{
    IPM_timer_getrusage *grut = (IPM_timer_getrusage *) ipmt;

    if (getrusage(RUSAGE_SELF, &(grut->start)))
        IPM_error_getrusage("IPM_timer_start()");
}

/* Do timeval calculations in microseconds as this avoids the accuracy loss */
/* that would occur when converting 1e-6 to binary if working in seconds. */
#define TV2US(timeval) (((timeval).tv_sec * 1e6) + (timeval).tv_usec)

static void
IPM_timer_stop_getrusage_process(IPM_timer *ipmt)
{
    IPM_timer_getrusage *grut = (IPM_timer_getrusage *) ipmt;

    if (getrusage(RUSAGE_SELF, &(grut->stop)))
        IPM_error_getrusage("IPM_timer_stop()");

    /* Add user and system time together. */
    grut->t_us += TV2US(grut->stop.ru_utime) - TV2US(grut->start.ru_utime);
    grut->t_us += TV2US(grut->stop.ru_stime) - TV2US(grut->start.ru_stime);
}

static void
IPM_timer_stop_getrusage_user(IPM_timer *ipmt)
{
    IPM_timer_getrusage *grut = (IPM_timer_getrusage *) ipmt;

    if (getrusage(RUSAGE_SELF, &(grut->stop)))
        IPM_error_getrusage("IPM_timer_stop()");

    /* Just add on user time. */
    grut->t_us += TV2US(grut->stop.ru_utime) - TV2US(grut->start.ru_utime);
}

static void
IPM_timer_stop_getrusage_system(IPM_timer *ipmt)
{
    IPM_timer_getrusage *grut = (IPM_timer_getrusage *) ipmt;

    if (getrusage(RUSAGE_SELF, &(grut->stop)))
        IPM_error_getrusage("IPM_timer_stop()");

    /* Just add on system time. */
    grut->t_us += TV2US(grut->stop.ru_stime) - TV2US(grut->start.ru_stime);
}

static double
IPM_timer_read_getrusage(const IPM_timer *ipmt)
{
    const IPM_timer_getrusage *grut = (const IPM_timer_getrusage *) ipmt;

    return grut->t_us * 1e-6;
}

static double
IPM_timer_get_time_getrusage_process(void)
{
    struct rusage r;

    if (getrusage(RUSAGE_SELF, &r))
        IPM_error_getrusage("IPM_timer_get_time()");

    return (TV2US(r.ru_utime) + TV2US(r.ru_stime)) * 1e-6;
}

static double
IPM_timer_get_time_getrusage_user(void)
{
    struct rusage r;

    if (getrusage(RUSAGE_SELF, &r))
        IPM_error_getrusage("IPM_timer_get_time()");

    return TV2US(r.ru_utime) * 1e-6;
}

static double
IPM_timer_get_time_getrusage_system(void)
{
    struct rusage r;

    if (getrusage(RUSAGE_SELF, &r))
        IPM_error_getrusage("IPM_timer_get_time()");

    return TV2US(r.ru_stime) * 1e-6;
}

void
IPM_init_getrusage(void)
{
    /* Parse timer type from environment. */
    char *tp = getenv("IPM_timer_type");

    IPM_timer_module = "getrusage";
    IPM_timer_type = tp;

    if (!tp || !strcmp(tp, "process"))
    {
        IPM_timer_type = "process";
        IPM_timer_stop_p = IPM_timer_stop_getrusage_process;
        IPM_timer_get_time_p = IPM_timer_get_time_getrusage_process;
    }
    else if (!strcmp(tp, "user"))
    {
        IPM_timer_type = "user";
        IPM_timer_stop_p = IPM_timer_stop_getrusage_user;
        IPM_timer_get_time_p = IPM_timer_get_time_getrusage_user;
    }
    else if (!strcmp(tp, "system"))
    {
        IPM_timer_type = "system";
        IPM_timer_stop_p = IPM_timer_stop_getrusage_system;
        IPM_timer_get_time_p = IPM_timer_get_time_getrusage_system;
    }
    else
        IPM_error("IPM_init(): No such IPM_timer_type.");

    /* Initialize resolution and range variables. */
    IPM_timer_resolution_val = IPM_GETRUSAGE_RESOLUTION;
    IPM_timer_range_val = IPM_GETRUSAGE_RANGE;

    /* Install remaining function pointers. */
    IPM_timer_clear_p = IPM_timer_clear_getrusage;
    IPM_timer_start_p = IPM_timer_start_getrusage;
    IPM_timer_read_p = IPM_timer_read_getrusage;
  
    /* Make sure opaque struct large enough. */
    assert(sizeof(IPM_timer_getrusage) <= IPM_TIMER_STRUCT_MAX_SIZE);
}
