const char* IPM_timer_t0vcount_rcsid = "$Header: /u/spert/src/IPM/RCS/IPM_timer_t0vcount.c,v 1.12 1996/08/07 23:37:12 krste Exp $";

/*
    Code to use the T0 on-chip vcount register on a SPERT board to
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

    Authors: Krste Asanovic <krste@cs.berkeley.edu>
             David Johnson <davidj@icsi.berkeley.edu>
*/

#include <assert.h>
#include <stdio.h>

#include "IPM_timer.h"
#include "IPM_timer_internal.h"
#include "IPM_timer_t0vcount.h"

/******************************************************************/
/* These parameters should really be retreived via a system call. */
/******************************************************************/

/* Assume 40 MHz boards for now (25 ns cycle time). */
#ifndef IPM_T0VCOUNT_RESOLUTION
#  define IPM_T0VCOUNT_RESOLUTION (25e-9)
#endif

typedef struct
{
    unsigned long start;
    unsigned long total;
} IPM_timer_t0vcount;

/*
    This routine reads the vcount register from coprocessor 2.  This
    holds a counter that is incremented once per clock cycle.
*/

#ifdef __GNUC__
static unsigned long
t0clock(void)
{
    register unsigned long res;

    __asm__ volatile ("cfvu %0, $1" : "=r" (res) : );
    return res;
}
#else
#error "Can only be compiled with gcc"
#endif

static void
IPM_timer_clear_t0vcount(IPM_timer *ipmt)
{
    IPM_timer_t0vcount *s0vc = (IPM_timer_t0vcount *) ipmt;

    s0vc->total = 0UL;
}

static void
IPM_timer_start_t0vcount(IPM_timer *ipmt)
{
    IPM_timer_t0vcount *s0vc = (IPM_timer_t0vcount *) ipmt;

    s0vc->start = t0clock();
}

static void
IPM_timer_stop_t0vcount(IPM_timer *ipmt)
{
    IPM_timer_t0vcount *s0vc = (IPM_timer_t0vcount *) ipmt;

    const unsigned long stop = t0clock(); /* No need to store this. */

    s0vc->total += (stop - s0vc->start);
}

static double
IPM_timer_read_t0vcount(const IPM_timer *ipmt)
{
    const IPM_timer_t0vcount *s0vc = (const IPM_timer_t0vcount *) ipmt;

    return IPM_timer_resolution_val * (double) s0vc->total;
}

static double
IPM_timer_get_time_t0vcount()
{
    return IPM_timer_resolution_val * (double) t0clock();
}

void
IPM_init_t0vcount(void)
{
    /* Make sure opaque struct large enough. */
    assert(sizeof(IPM_timer_t0vcount) <= IPM_TIMER_STRUCT_MAX_SIZE);

    IPM_timer_module = "t0vcount";
    IPM_timer_type = "real";

    IPM_timer_resolution_val = IPM_T0VCOUNT_RESOLUTION;
    IPM_timer_range_val = IPM_T0VCOUNT_RESOLUTION * (double) 0xffffffffUL;

    /* Install function pointers. */
    IPM_timer_clear_p = IPM_timer_clear_t0vcount;
    IPM_timer_start_p = IPM_timer_start_t0vcount;
    IPM_timer_stop_p = IPM_timer_stop_t0vcount;
    IPM_timer_read_p = IPM_timer_read_t0vcount;
    IPM_timer_get_time_p = IPM_timer_get_time_t0vcount;
}
