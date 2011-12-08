/*
    $Header: /u/spert/src/IPM/RCS/IPM_timer_gettimer.h,v 1.2 1995/10/14 23:13:06 krste Exp $

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

#ifndef IPM_timer_gettimer_h_
#define IPM_timer_gettimer_h_

void IPM_init_gettimer(void);

#endif
