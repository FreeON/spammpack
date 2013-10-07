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
#include "utilities.h"
#include "index.h"

#include <getopt.h>

/** The main program.
 *
 * @param msg The command line arguments.
 */
BCSRInfo::BCSRInfo (CkArgMsg *msg)
{
  char *MM_filename = NULL;
  int c;
  const char *short_options = "hw:";
  const struct option long_options[] = {
    { "help",     no_argument,        NULL, 'h' },
    { "write-MM", required_argument,  NULL, 'w' },
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
        CkPrintf("\n");
        CkPrintf("{ -h | --help }             This help\n");
        CkPrintf("{ -w | --write-MM } FILE    Write matrix in MatrixMarket format to FILE\n");
        CkExit();
        break;

      case 'w':
        MM_filename = strdup(optarg);
        break;

      default:
        ABORT("unknown command line option\n");
        break;
    }
  }

  if(optind < msg->argc)
  {
    BCSR A(msg->argv[optind]);
    double F_min, F_max;
    A.getSpectralBounds(0, &F_min, &F_max);
    printf("spectral bounds Gershgorin: [ %e, %e ]\n", F_min, F_max);
    A.getSpectralBounds(1, &F_min, &F_max);
    printf("spectral bounds eigensolve: [ %e, %e ]\n", F_min, F_max);
    A.toStr();

    if(MM_filename != NULL)
    {
      A.toMM(MM_filename);
    }

    CkExit();
  }

  else
  {
    ABORT("missing filename\n");
  }
}

#include "bcsrinfo.def.h"
