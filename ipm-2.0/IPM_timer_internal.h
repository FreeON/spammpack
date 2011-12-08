/*
    $Header: /u/spert/src/IPM/RCS/IPM_timer_internal.h,v 1.4 1996/08/08 00:43:44 krste Exp $

    Internal header file for IPM_timer.

    IPM timer package provides portable, high resolution interval timing.
*/

/*
    "Copyright (c) 1996 The Regents of the University of California.  
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

#ifndef IPM_timer_internal_h_
#define IPM_timer_internal_h_

extern const char* IPM_version;
extern const char* IPM_version_str;

extern double IPM_timer_resolution_val;
extern double IPM_timer_range_val;

extern char* IPM_timer_module;
extern char* IPM_timer_type;

extern void IPM_error(const char*);
extern void IPM_system_error(const char*);

#endif
