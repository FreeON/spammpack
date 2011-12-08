const char* IPM_timer_syssgi_rcsid = "$Header: /u/spert/src/IPM/RCS/IPM_timer_syssgi.c,v 1.18 1996/08/08 00:43:46 krste Exp $";

/*
    Code to use the hardware real time counter on SGI systems.

    This release should work across nearly all SGI platforms, thanks to
    example code from Dave Cortesi <cortesi@barchester.wpd.sgi.com> at SGI.

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

/*
    Enhanced with code extracted from the REACT/Pro Programmer's Guide,
    part #007-2499-001,
    Copyright 1994 Silicon Graphics, Inc.

    Supplied by Dave Cortesi.
*/

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/syssgi.h>
#include <sys/mman.h>
#include <unistd.h>         /* for getpagesize() */
#include <sgidefs.h>        /* for __psint_t, __uint32_t, and __int64_t */
#include <invent.h>         /* getinvent() and related constants */
#include <stdio.h>

#include "IPM_timer.h"
#include "IPM_timer_internal.h"
#include "IPM_timer_syssgi.h"

typedef struct
{
    union {
        __uint32_t ui32;
        __uint64_t ui64;
    } start;
    union {
        __uint32_t ui32;
        __uint64_t ui64;
    } total;
} IPM_timer_syssgi;

static volatile void* clock_p;

static int timer_bits;          /* Number of timer bits. */

/*

  The *_32 calls are used when timer is 24 or 32 bits.
  The *_64 calls are used when timer is 64 bits.

*/

static void
IPM_timer_clear_syssgi_32(IPM_timer *ipmt)
{
    IPM_timer_syssgi *sst = (IPM_timer_syssgi *) ipmt;

    sst->total.ui32 = 0;
}

static void
IPM_timer_start_syssgi_32(IPM_timer *ipmt)
{
    IPM_timer_syssgi *sst = (IPM_timer_syssgi *) ipmt;

    sst->start.ui32 = *((__uint32_t*) clock_p);
}

static void
IPM_timer_stop_syssgi_32(IPM_timer *ipmt)
{
    IPM_timer_syssgi *sst = (IPM_timer_syssgi *) ipmt;

    const __uint32_t stop = *((__uint32_t*) clock_p);

    sst->total.ui32 += (stop - sst->start.ui32);
}

static double
IPM_timer_read_syssgi_32(const IPM_timer *ipmt)
{
    const IPM_timer_syssgi *sst = (const IPM_timer_syssgi *) ipmt;

    return IPM_timer_resolution_val * (double) sst->total.ui32;
}

static double
IPM_timer_get_time_syssgi_32(void)
{
    return IPM_timer_resolution_val * (double)(*((__uint32_t*) clock_p));
}

static void
IPM_timer_clear_syssgi_64(IPM_timer *ipmt)
{
    IPM_timer_syssgi *sst = (IPM_timer_syssgi *) ipmt;

    sst->total.ui64 = 0;
}

static __uint64_t
read_clock_64(void)
{
#if (_MIPS_SZPTR == 64) /* 64-bit program, 64-bit timer */
            return *((__uint64_t*) clock_p);
#else /* 32-bit program, 64-bit timer */
            __uint32_t msw_start = ((__uint32_t*)clock_p)[0];
            __uint32_t lsw = ((__uint32_t*)clock_p)[1];

            __uint32_t msw_end = ((__uint32_t*)clock_p)[0];

            while (msw_end != msw_start) /* Wrap-around during read. */
            {
                msw_start = msw_end;
                lsw = ((__uint32_t*)clock_p)[1];
                msw_end =  ((__uint32_t*)clock_p)[0];
            }

            /* lsw read while msw didn't change. */

            return (((__uint64_t) msw_start)<<32)|((__uint64_t)lsw);
#endif
}

static void
IPM_timer_start_syssgi_64(IPM_timer *ipmt)
{
    IPM_timer_syssgi *sst = (IPM_timer_syssgi *) ipmt;

    sst->start.ui64 = read_clock_64();
}

static void
IPM_timer_stop_syssgi_64(IPM_timer *ipmt)
{
    IPM_timer_syssgi *sst = (IPM_timer_syssgi *) ipmt;

    const __uint64_t stop = read_clock_64();

    sst->total.ui64 += (stop - sst->start.ui64);
}

static double
IPM_timer_read_syssgi_64(const IPM_timer *ipmt)
{
    const IPM_timer_syssgi *sst = (const IPM_timer_syssgi *) ipmt;

    return IPM_timer_resolution_val * (double) sst->total.ui64;
}

static double
IPM_timer_get_time_syssgi_64(void)
{
    return IPM_timer_resolution_val * (double) read_clock_64();
}

/*
|| getIPNumber()
||
|| This function returns the IPnn number -- the CPU board model number --
|| of the current system.  For example in an Indy it returns 22 (IP22),
|| while in a Challenge-XL it returns 19 (IP19).  Used to determine number
|| of bits in the available timer.
*/

static unsigned short
getIPNumber()
{
    inventory_t *pitem;
    unsigned short ipno = 0; /* default meaning ?? */

    setinvent();
    while( (pitem = getinvent()) )
    {
        if ((pitem->inv_class == INV_PROCESSOR)
            &&  (pitem->inv_type  == INV_CPUBOARD) )
        {
            switch (pitem->inv_state)
            {
            case INV_IP4BOARD: ipno = 4; break;
            case INV_IP5BOARD: ipno = 5; break;
            case INV_IP6BOARD: ipno = 6; break;
            case INV_IP7BOARD: ipno = 7; break;
            case INV_IP9BOARD: ipno = 9; break;
            case INV_IP12BOARD: ipno = 12; break;
            case INV_IP17BOARD: ipno = 17; break;
            case INV_IP20BOARD: ipno = 20; break;
            case INV_IP19BOARD: ipno = 19; break;
            case INV_IP22BOARD: ipno = 22; break;
            case INV_IP21BOARD: ipno = 21; break;
            case INV_IP26BOARD: ipno = 26; break;
            default: break;
            }
            break;
        }
    }
    endinvent();
    return ipno;
}

static void
IPM_init_counter(void)
{
    __psint_t paddr;            /* Return value from syssgi(). */
    int mem_fd;                 /* File descriptor of system memory. */
    int pagesize;               /* Size of system pages. */
    __psint_t page_addr;        /* Page aligned address of counter. */
    __psint_t clock_addr_lo;    /* Low bits of clock address. */
    void *vaddr;                /* Virtual page address returned by mmap(). */
    __uint32_t picoseconds;     /* Resolution in picoseconds. */
    __int64_t max_val;          /* Maximum value in timer register. */

    /* Check for cycle counter. */
    paddr = syssgi(SGI_QUERY_CYCLECNTR, &picoseconds);
    if (paddr == (__psint_t) -1)
        IPM_system_error("IPM_init(): No hardware cycle counter present.");

    /* fprintf(stderr, "paddr=0x%lx\n", (unsigned long) paddr); */

    mem_fd = open("/dev/mmem", O_RDONLY);
    if (mem_fd == -1)
        IPM_system_error("IPM_init(): open()");

    /* Map cycle counter into page-aligned user memory. */
    pagesize = getpagesize();
    clock_addr_lo = paddr & (pagesize - 1);
    page_addr = paddr - clock_addr_lo;
    vaddr = mmap(0,             /* Don't care where in virtual space. */
                 pagesize,      /* Only need one page. */
                 PROT_READ,     /* Only need to read. */
                 MAP_SHARED,    /* Simplest option. */
                 mem_fd,        /* From /dev/mem open command. */
                 (off_t) page_addr); /* Clock's physical address. */

    if ((__psint_t) vaddr == (__psint_t) -1)
        IPM_system_error("IPM_init(): mmap():");

    /* fprintf(stderr, "vaddr=0x%lx\n", (unsigned long) vaddr); */

    clock_p = (void *) ((__psint_t) vaddr + clock_addr_lo);

    /* Convert picoseconds to seconds. */
    IPM_timer_resolution_val = picoseconds * 1e-12;

    switch (getIPNumber())
    {
    case 12: /* IP12 (Indigo R3K) */
        timer_bits = 24;
        max_val = (__uint64_t) 0x00ffffff;
        break;
    case 19: /* IP19 (Challenge, Onyx) */
    case 21:
        timer_bits = 64;
        max_val = ~((__uint64_t)0);
        break;
    default: /* everything else currently is... */
        timer_bits = 32;
        max_val = (__uint64_t) 0xffffffff; /* ...32 bits, Dave Olson says */
        break;
    }

    IPM_timer_range_val = (double) max_val * IPM_timer_resolution_val;
}

void
IPM_init_syssgi(void)
{
    /* Make sure opaque struct large enough. */
    assert(sizeof(IPM_timer_syssgi) <= IPM_TIMER_STRUCT_MAX_SIZE);

    IPM_timer_module="syssgi";
    IPM_timer_type="real";

    /* Set up counter. */
    IPM_init_counter();

    /* Install function pointers. */
    if (timer_bits<=32)
    {
        IPM_timer_clear_p = IPM_timer_clear_syssgi_32;
        IPM_timer_start_p = IPM_timer_start_syssgi_32;
        IPM_timer_stop_p = IPM_timer_stop_syssgi_32;
        IPM_timer_read_p = IPM_timer_read_syssgi_32;
        IPM_timer_get_time_p = IPM_timer_get_time_syssgi_32;
    }
    else
    {
        IPM_timer_clear_p = IPM_timer_clear_syssgi_64;
        IPM_timer_start_p = IPM_timer_start_syssgi_64;
        IPM_timer_stop_p = IPM_timer_stop_syssgi_64;
        IPM_timer_read_p = IPM_timer_read_syssgi_64;
        IPM_timer_get_time_p = IPM_timer_get_time_syssgi_64;
    }
}
