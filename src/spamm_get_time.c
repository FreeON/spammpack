/** @file */

/*
     This code is part of the FreeON suite of programs for linear scaling
     electronic structure theory and ab initio molecular dynamics.

     Copyright (2004). The Regents of the University of California. This
     material was produced under U.S. Government contract W-7405-ENG-36
     for Los Alamos National Laboratory, which is operated by the University
     of California for the U.S. Department of Energy. The U.S. Government has
     rights to use, reproduce, and distribute this software.  NEITHER THE
     GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
     OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by the
     Free Software Foundation; either version 3 of the License, or (at your
     option) any later version. Accordingly, this program is distributed in
     the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. See the GNU General Public License at www.gnu.org for details.

     While you may do as you like with this software, the GNU license requires
     that you clearly mark derivative software.  In addition, you are encouraged
     to return derivative works to the FreeON group for review, and possible
     dissemination in future releases.
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>

/** Get the current time.
 */
double
spamm_get_time ()
{
#ifdef TIMER_GETRUSAGE
  struct rusage now;

  if (getrusage(RUSAGE_SELF, &now) != 0)
  {
    printf("error running getrusage()\n");
    exit(1);
  }
  return (now.ru_utime.tv_sec*1e6+now.ru_utime.tv_usec)/1.0e6;
#elif defined(TIMER_GETTIMEOFDAY)
  struct timeval now;

  if (gettimeofday(&now, NULL) != 0)
  {
    printf("error running gettimeofday()\n");
    exit(1);
  }
  return (now.tv_sec*1e6+now.tv_usec)/1.0e6;
#else
  printf("no timer configured\n");
  exit(1);
#endif
}

/** @private
 *
 * Underscore version of spamm_get_time().
 */
double
spamm_get_time_ ()
{
  return spamm_get_time();
}

/** @private
 *
 * Underscore version of spamm_get_time().
 */
double
spamm_get_time__ ()
{
  return spamm_get_time();
}
