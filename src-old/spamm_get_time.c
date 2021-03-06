/** @file
 *
 * @copyright
 *
 * Copyright (c) 2015, Los Alamos National Laboratory
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @author Matt Challacombe matt.challacombe@freeon.org
 * @author Nicolas Bock nicolasbock@freeon.org
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>

/** Get the current time.
 *
 * @param timer The time.
 */
void
spamm_get_time (double *timer)
{
#ifdef TIMER_GETRUSAGE
  struct rusage now;

  if (getrusage(RUSAGE_SELF, &now) != 0)
  {
    printf("error running getrusage()\n");
    exit(1);
  }
  *timer = (now.ru_utime.tv_sec*1e6+now.ru_utime.tv_usec)/1.0e6;
#elif defined(TIMER_GETTIMEOFDAY)
  struct timeval now;

  if (gettimeofday(&now, NULL) != 0)
  {
    printf("error running gettimeofday()\n");
    exit(1);
  }
  *timer = (now.tv_sec*1.0e6+now.tv_usec)/1.0e6;
#else
  printf("no timer configured\n");
  exit(1);
#endif
}
