const char* IPM_timer_gettimer_rcsid = "$Header: /u/spert/src/IPM/RCS/IPM_timer_gettimer.c,v 1.15 1996/08/08 00:43:43 krste Exp $";

/*
    Code that uses AIX gettimer() to implement an IPM timer.

    Based on suggestion in IBM's:
    Optimization and Tuning Guide for Fortran, C, and C++ (order # SC09-1705)
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
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

#include "IPM_timer.h"
#include "IPM_timer_internal.h"
#include "IPM_timer_gettimer.h"

typedef struct
{
    struct timestruc_t start;   /* Start values. */
    struct timestruc_t stop;    /* Stop values. */
    double t_ns;                /* Accumulated nanoseconds, stored in */
                                /* double. */
} IPM_timer_gettimer;

static void
IPM_timer_clear_gettimer(IPM_timer *ipmt)
{
    IPM_timer_gettimer *gt = (IPM_timer_gettimer *) ipmt;

    gt->t_ns = 0.0;              /* Clear accumulated time. */
}

static void
IPM_timer_start_gettimer(IPM_timer *ipmt)
{
    IPM_timer_gettimer *gt = (IPM_timer_gettimer *) ipmt;

    if (gettimer(TIMEOFDAY, &(gt->start)))
        IPM_system_error("IPM_timer_start(): gettimer()");
}

/* Do timestruc_t calculations in nanoseconds as this avoids the accuracy */
/* loss that would occur when converting 1e-9 to binary if working in    */
/* seconds. */
#define TS2NS(ts) \
    (((ts).tv_sec * 1e9) + (ts).tv_nsec)

static void
IPM_timer_stop_gettimer(IPM_timer *ipmt)
{
    IPM_timer_gettimer *gt = (IPM_timer_gettimer *) ipmt;

    if (gettimer(TIMEOFDAY, &(gt->stop)))
        IPM_system_error("IPM_timer_stop(): gettimer()");
    gt->t_ns += TS2NS(gt->stop) - TS2NS(gt->start);
}

static double
IPM_timer_read_gettimer(const IPM_timer *ipmt)
{
    const IPM_timer_gettimer *gt = (const IPM_timer_gettimer *) ipmt;
    
    return gt->t_ns * 1e-9;
}

static double
IPM_timer_get_time_gettimer(void)
{
    struct timestruc_t t;
    
    if (gettimer(TIMEOFDAY, &t))
        IPM_system_error("IPM_timer_get_time(): gettimer()");

    return TS2NS(t) * 1e-9;
}

void
IPM_init_gettimer(void)
{
    struct timestruc_t res;
    struct timestruc_t max_value;

    /* Make sure opaque struct large enough. */
    assert(sizeof(IPM_timer_gettimer) <= IPM_TIMER_STRUCT_MAX_SIZE);

    IPM_timer_module = "gettimer";
    IPM_timer_type = "real";

    if (restimer(TIMEOFDAY, &res, &max_value))
        IPM_system_error("IPM_init(): restimer()");

    IPM_timer_resolution_val = TS2NS(res) * 1e-9;
    IPM_timer_range_val = TS2NS(max_value) * 1e-9;

    /* Install function pointers. */
    IPM_timer_clear_p = IPM_timer_clear_gettimer;
    IPM_timer_start_p = IPM_timer_start_gettimer;
    IPM_timer_stop_p = IPM_timer_stop_gettimer;
    IPM_timer_read_p = IPM_timer_read_gettimer;
    IPM_timer_get_time_p = IPM_timer_get_time_gettimer;
}
