const char* IPM_timer_stats_rcsid = "$Header: /u/spert/src/IPM/RCS/IPM_timer_stats.c,v 1.11 1996/08/06 06:15:17 krste Exp $";

/*
    Measures properties of a timer, using the timer itself.

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
*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifndef EXIT_SUCCESS
# define EXIT_SUCCESS (0)
#endif
#ifndef EXIT_FAILURE
# define EXIT_FAILURE (0)
#endif

#include "IPM_timer.h"

/*
  Stats routines provide simple statistical measures on vector of doubles.
  Includes routine to find minimum non-zero difference between any two
  values.
*/

typedef struct {
    double min;                 /* Minimum value. */
    double percentile10;        /* 10-percentile. */
    double quartile1;           /* 25-percentile (first quartile). */
    double median;              /* Median value (second quartile). */
    double quartile3;           /* 75-percentile (third quartile). */
    double percentile90;        /* 90-percentile. */
    double max;                 /* Maximum value. */
    double siqr;                /* Semi-InterQuartile Range. */

    double mean;                /* Mean value. */
    double stddev;              /* Standard Deviation. */
    double cov;                 /* Coefficient of Variance. */

    double min_delta;           /* Minimum non-zero difference between
                                   samples. */
} Stats;

static int
double_compar(const void *ap, const void *bp)
{
    const double a = *((double *) ap);
    const double b = *((double *) bp);

    if (a < b)
        return -1;
    else if (a > b)
        return +1;
    else
        return 0;
}

static void
Stats_calc(long n, double *results, Stats* s)
{
    long i;
    double sum;
    double sum2;
    double mean;

    double min_delta;

    assert(n>=2);

    /* Sort results to get quantiles and min and max, and to simplify
       min_delta calculation. */
    qsort((void *) results, n, sizeof(double), double_compar);

    s->min = results[0];
    s->percentile10 = results[n/10];
    s->quartile1 = results[n/4];
    s->median = results[n/2];
    s->quartile3 = results[3*n/4];
    s->percentile90 = results[9*n/10];
    s->max = results[n-1];

    /* Get Semi-Interquartile Range. */
    s->siqr = (s->quartile3 - s->quartile1)/2;

    /* Calculate mean and stddev in two passes to avoid precision hassles. */
    sum = 0.0;
    for (i=0; i<n; i++)
    {
        sum += results[i];
    }
    mean = sum/n;               /* Keep mean in local to speed next loop. */

    sum2 = 0.0;
    for (i=0; i<n; i++)
    {
        double diff = results[i] - mean;
        sum2 += diff * diff;
    }

    s->mean = mean;
    s->stddev = sqrt(sum2/(n - 1));
    s->cov = s->stddev/s->mean;

    /* Find smallest non-zero difference between two readings. */
    min_delta = 0.0;
    for (i=1; i<n; i++)
    {
        const double delta = results[i] - results[i-1];

        if (delta > 0.0)
        {
            if (min_delta > 0.0)
            {
                if (delta < min_delta)
                    min_delta = delta;
            }
            else
                min_delta = delta;
        }
    }

    s->min_delta = min_delta;
}

static void
Stats_scale(double scale, Stats *stats)
{
    stats->min *= scale;
    stats->percentile10 *= scale;
    stats->quartile1 *= scale;
    stats->median *= scale;
    stats->quartile3 *= scale;
    stats->percentile90 *= scale;
    stats->max *= scale;
    stats->siqr *= scale;
    
    stats->mean *= scale;
    stats->stddev *= scale;
    /* Don't scale cov. */
    /* Don't scale min_delta. */
}

static void
Stats_print(Stats* stats)
{
    (void) printf("Minimum time  %4e s.\n", stats->min);
    (void) printf("10-percentile %4e s.\n", stats->percentile10);
    (void) printf("25-percentile %4e s.\n", stats->quartile1);
    (void) printf("Median time   %4e s.\n", stats->median);
    (void) printf("75-percentile %4e s.\n", stats->quartile3);
    (void) printf("90-percentile %4e s.\n", stats->percentile90);
    (void) printf("Maximum time  %4e s.\n", stats->max);
    (void) printf("SIQR          %4e s.\n", stats->siqr);
    (void) printf("\n");
    (void) printf("Mean time     %4e s.\n", stats->mean);
    (void) printf("Stddev.       %4e s.\n", stats->stddev);
    (void) printf("COV           %4e.\n",   stats->cov);
    (void) printf("\n");

    if (stats->min_delta > 0.0)
    {
        (void) printf(
            "Smallest delta measured was        %4e s.\n",
            stats->min_delta);
    }
    else
    {
        (void) printf("No delta could be measured.\n");
    }
    (void) printf(
        "Value of IPM_timer_resolution() is %4e s.\n", IPM_timer_resolution());

}

/*
  Try to measure overhead experienced by one IPM_timer_start/stop pair
  using just the timer itself.
*/

static void
one_shot_measure(long n, double *results)
{
    IPM_timer t;
    long i;

    for (i=0; i<n; i++)
    {
        IPM_timer_clear(&t);
        IPM_timer_start(&t);
        IPM_timer_stop(&t);
        results[i] = IPM_timer_read(&t);
    }
}

static void
one_shot_run(long n)
{
    double *results;
    Stats stats;

    (void) printf("----------------------------------------------------------------------------\n");
    (void) printf("Measuring overhead for a single "
                  "IPM_timer_start/stop pair.\n");
    (void) printf("Can be inaccurate when timer increments large "
                  "compared to overhead.\n");

    results = (double *) malloc(n * sizeof(double));
    if (results==NULL)
    {
        (void) fprintf(stderr, "Couldn't allocate results buffer.\n");
        exit(EXIT_FAILURE);
    }

    (void) printf("Taking measurements.\n");
    one_shot_measure(n, results);
    (void) printf("Finished taking measurements.\n");

    (void) printf("Calculating statistics over %ld samples.\n", n);
    Stats_calc(n, results, &stats);
    Stats_print(&stats);
    (void) printf("----------------------------------------------------------------------------\n");

    free(results);
}

/*
  Following routines support measurement of total overhead of one
  IPM_timer_start call and one IPM_timer_stop call.

  This is different than the overhead for an IPM_timer measurement.  Where
  the system timer call dominates execution time, the total overhead is
  about twice the overhead experienced by an IPM_timer measurement.
*/

static void
total_overhead_measure(long n_outer, long n_inner, double *results)
{
    IPM_timer t;
    IPM_timer dummy_timer;
    long outer;
    long inner;

    for (outer=0; outer<n_outer; outer++)
    {
        IPM_timer_clear(&t);
        IPM_timer_clear(&dummy_timer);
        IPM_timer_start(&t);
        for (inner=0; inner<n_inner; inner++)
        {
            IPM_timer_start(&dummy_timer);
            IPM_timer_stop(&dummy_timer);
        }
        IPM_timer_stop(&t);

        results[outer] = IPM_timer_read(&t);
    }
}

#define N_OUTER_DEFAULT (100)

#define MIN_INNER_TIME (0.01)
/* 10 outer loops used during inner loop calibration. */
#define CALIBRATE_OUTER (10)
/* Minimum multiple of timer resolution. */
#define CALIBRATE_RES_MULT (100)

static long
total_overhead_find_n_inner()
{
    /*
       Finds a good value for n_inner.  Run inner loop till run time is at
       least MIN_INNER_TIME, or at least CALIBRATE_RES_MULT times
       resolution, whichever is largest.
       */

    const double resolution = IPM_timer_resolution();
    const double min_res_time = resolution * CALIBRATE_RES_MULT;
    const double min_time = min_res_time > MIN_INNER_TIME
        ? min_res_time : MIN_INNER_TIME;

    double results[CALIBRATE_OUTER];

    long n_inner;

    (void) printf("Running auto-calibrate to pick n_inner.\n");

    n_inner = 10;
    while (1)
    {
        long i;
        double min;

        (void) printf("n_inner = %10ld,", n_inner);
        total_overhead_measure(CALIBRATE_OUTER, n_inner, results);

        /* Find minimum in results array. */
        min = results[CALIBRATE_OUTER - 1]; /* Last one probably smaller. */
        for (i = 0; i < CALIBRATE_OUTER - 1; i++)
        {
            const double r = results[i];
            if (r < min)
                min = r;
        }

        (void) printf(" minimum inner loop time = %4e s.\n", min);

        
        if (min > min_time)
        {
            n_inner = (n_inner * min_time/min) + 10;
            break;
        }
        else if (min > 0.0)
        {
            /* Try to guess size needed, add 10% to help ensure big enough. */
            n_inner = (n_inner * min_time/min) * 1.10;
        }
        else
            n_inner *= 10;      /* Go up 10 times, if min==0.0 */
    }

    (void) printf("Auto-calibrate selecting n_inner = %ld.\n", n_inner);

    return n_inner;
}

static void
total_overhead_run(long n_outer, /* Number of measurements to take. */
                   long n_inner /* Number inside inner loop, if zero run
                                   automatic routine to determine good
                                   value. */
    )
{
    double *results;
    Stats stats;

    (void) printf("----------------------------------------------------------------------------\n");
    (void) printf("Measuring total overhead of "
                  "one IPM_timer_start and one IPM_timer_stop call.\n");
    (void) printf("This is approximately twice the overhead seen by a "
                  "single IPM_timer measurement.\n");

    results = (double *) malloc(n_outer * sizeof(double));
    if (results==NULL)
    {
        (void) fprintf(stderr, "Couldn't allocate results buffer.\n");
        exit(EXIT_FAILURE);
    }

    if (!n_inner)
    {
        n_inner = total_overhead_find_n_inner();
    }

    (void) printf("Starting measurements.\n");
    total_overhead_measure(n_outer, n_inner, results);
    (void) printf("Finished measurements.\n");

    (void) printf("Calculating statistics over %ld samples, "
                  "each averaged over %ld calls.\n", n_outer, n_inner);
    Stats_calc(n_outer, results, &stats);
    Stats_scale(1.0/n_inner, &stats);
    Stats_print(&stats);
    (void) printf("----------------------------------------------------------------------------\n");

    free(results);
}

static void
usage(int exit_status)
{
    (void) fprintf(stderr,
"IPM_timer_stats [-help] [-n <#outer>] [-i #inner]\n"
"    -help       Show this message.\n"
"    -n <#outer> Use this many outer loop timer readings.  Default is 100.\n"
"    -i <#inner> Call routine this many times in inner loop.  By default,\n"
"                the program automatically determines this parameter.\n");

    exit(exit_status);
}

static void
parse_long(const char*const s, long *const lp)
{
    size_t len = strlen(s);
    char *ptr;

    *lp = strtol(s, &ptr, 0);

    if (ptr != (s+len))
    {
        (void) fprintf(stderr, "ERROR: Value %s not an integer.\n", s);
        usage(EXIT_FAILURE);
    }
}

int
main(int argc, char *argv[])
{
    long n_outer = N_OUTER_DEFAULT;
    long n_inner = 0;           /* Zero means find n_inner automatically. */

    ++argv;                     /* Move past program name. */

    while (--argc)
    {
        if (!strcmp(*argv, "-n"))
        {
            if (--argc)
            {
                parse_long(*(++argv), &n_outer);
                ++argv;
            }
            else
            {
                (void) fprintf(stderr,
                               "ERROR: No outer loop count specified.\n");
                usage(EXIT_FAILURE);
            }
        }
        else if (!strcmp(*argv, "-i"))
        {
            if (--argc)
            {
                parse_long(*(++argv), &n_inner);
                ++argv;
            }
            else
            {
                (void) fprintf(stderr,
                               "ERROR: No inner loop count specified.\n");
                usage(EXIT_FAILURE);
            }
        }
        else if (!strcmp(*argv, "-help"))
        {
            usage(EXIT_SUCCESS);
        }
        else
        {
            (void) fprintf(stderr, "ERROR: Unrecognized option: %s.\n", *argv);
            usage(EXIT_FAILURE);
        }
    }

    one_shot_run(n_outer);
    total_overhead_run(n_outer, n_inner);

    return EXIT_SUCCESS;
}
