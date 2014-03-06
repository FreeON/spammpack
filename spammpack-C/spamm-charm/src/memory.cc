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

/** Return the current virtual memory usage of the caller.
 *
 * @return The current usage in kiB.
 */
int Memory::get_virtual (void)
{
  int result = 0;

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
    char *found_it = strstr(line, "VmSize");
    if(found_it != NULL)
    {
      if(sscanf(line, "VmSize: %d kB", &result) != 1)
      {
        printf("can not parse VmSize\n");
        exit(1);
      }
      break;
    }
  }
  fclose(stream);
  return result;
}
