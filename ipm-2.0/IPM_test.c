const char *IPM_test_rcsid = "$Header: /u/spert/src/IPM/RCS/IPM_test.c,v 1.16 1996/08/01 23:16:27 krste Exp $";

/*
    Test wrapper for IPM timer routines.
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

#include <stdio.h>

#include "IPM_timer.h"

static char foo[100];

/* Timing loop parameters. */
#define CLEAR_LOOPS (2)
#define ACC_LOOPS (10)
#define INNER_LOOPS (10000)

int
main()
{
    IPM_timer t;
    long k;
    long j;
    long i;
    long sum = 0;
    volatile char* p = foo;
    char buf[IPM_TIMER_REPORT_LEN+20];

    /* Check out range and resolution calls. */
    printf("Timer range %3g s, resolution %4g s\n",
           IPM_timer_range(),
           IPM_timer_resolution());

    /* Test IPM_timer_format(). */
    printf("IPM_timer_format returns %s.\n", IPM_timer_format());

    /* Try interval timing and timer report. */
    for (k=0; k<CLEAR_LOOPS; k++)
    {
        IPM_timer_clear(&t);

        for (j=0; j<ACC_LOOPS; j++)
        {
            IPM_timer_start(&t);

            for (i=0; i<INNER_LOOPS; i++)
            {
                sum += p[(i % 12) + 2];
            }

            IPM_timer_stop(&t);

            IPM_timer_report("testing", &t, buf);

            printf("%2ld: Timer reports\n%s", j, buf);
        }
        printf("\n");
    }

    /* Try out IPM_timer_read call. */
    for (k=0; k<CLEAR_LOOPS; k++)
    {
        IPM_timer_clear(&t);

        for (j=0; j<ACC_LOOPS; j++)
        {
            IPM_timer_start(&t);

            for (i=0; i<INNER_LOOPS; i++)
            {
                sum += p[(i % 12) + 2];
            }

            IPM_timer_stop(&t);

            printf("%2ld: IPM_timer_read = %g s.\n", j, IPM_timer_read(&t));
        }
        printf("\n");
    }

    /* Try out IPM_timer_get_time call. */

    for (k=0; k<CLEAR_LOOPS; k++)
    {
        double t;
        double last_t = 0;

        for (j=0; j<ACC_LOOPS; j++)
        {
            t = IPM_timer_get_time();
            printf("%2ld: IPM_timer_get_time = %20.18g s, diff = %20.18g s.\n",
                   j,
                   t,
                   t - last_t);

            last_t = t;

            for (i=0; i<INNER_LOOPS; i++)
            {
                sum += p[(i % 12) + 2];
            }
        }

        printf("\n");
    }

    return 0;
}
