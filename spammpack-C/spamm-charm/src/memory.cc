/** @file
 *
 * Implementation of the Memory class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "memory.h"

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

/** Parse the process' status in the /proc FS.
 *
 * @param VmSize The current virtual memory size.
 * @param VmPeak The peak virtual memory usage.
 */
void Memory::parse_proc (int *VmSize, int *VmPeak)
{
  /* Get PID of caller (well it's really the PID of this process here). */
  pid_t my_pid = getpid();
  char filename[2000];

  /* Open the proper file in the /proc FS. */
  snprintf(filename, 2000, "/proc/%d/status", my_pid);
  int fd = open(filename, O_RDONLY);

  /* Read the file. */
  FILE *stream = fdopen(fd, "r");
  char line[2000];
  while(fgets(line, 2000, stream) != NULL)
  {
    if(VmSize != NULL)
    {
      char *found_it = strstr(line, "VmSize");
      if(found_it != NULL)
      {
        if(sscanf(line, "VmSize: %d kB", VmSize) != 1)
        {
          printf("can not parse VmSize\n");
          exit(1);
        }
      }
    }

    if(VmPeak != NULL)
    {
      char *found_it = strstr(line, "VmPeak");
      if(found_it != NULL)
      {
        if(sscanf(line, "VmPeak: %d kB", VmPeak) != 1)
        {
          printf("can not parse VmPeak\n");
          exit(1);
        }
      }
    }
  }
  fclose(stream);
}

/** Return the current virtual memory usage of the caller.
 *
 * @return The current usage in kiB.
 */
int Memory::get_virtual (void)
{
  int VmSize;
  Memory::parse_proc(&VmSize, NULL);
  return VmSize;
}

/** Return the peak virtual memory usage of the caller.
 *
 * @return The peak usage in kiB.
 */
int Memory::get_peak_virtual (void)
{
  int VmPeak;
  Memory::parse_proc(NULL, &VmPeak);
  return VmPeak;
}
