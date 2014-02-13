#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

int
main (int argc, char **argv)
{
  char *F_filename = NULL;

  int c;
  const char short_options[] = "h";
  const struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
  };

  while((c = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
  {
    switch(c)
    {
      default:
        printf("illegal command line argument\n");
        exit(-1);
        break;
    };
  }

  if(F_filename == NULL)
  {
    printf("missing Fockian matrix file\n");
    exit(-1);
  }

  return 0;
}
