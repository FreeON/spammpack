/** @file
 *
 * A program to print out some information on a BCSR file.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "bcsrinfo.h"
#include "bcsr.h"
#include "logger.h"

#include <getopt.h>

BCSRInfo::BCSRInfo (CkArgMsg *msg)
{
  int c;
  const char *short_options = "h";
  const struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
  };

  initializeLogger();

  while((c = getopt_long(msg->argc, msg->argv, short_options, long_options,
          NULL)) != -1)
  {
    switch(c)
    {
      case 'h':
        CkPrintf("Usage:\n");
        CkExit();
        break;

      default:
        ABORT("unknown command line option\n");
        break;
    }
  }

  if(optind < msg->argc)
  {
    BCSR A(msg->argv[optind]);
    A.toStr();
    CkExit();
  }

  else
  {
    ABORT("missing filename\n");
  }
}

#include "bcsrinfo.def.h"
