#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

int
main (int argc, char **argv)
{
  const char *short_options = "c:b:h";
  const struct option long_options[] = {
    { "N_chunk", required_argument, NULL, 'c' },
    { "N_basic", required_argument, NULL, 'b' },
    { "help", no_argument, NULL, 'h' }
  };
  char c;
  int N_chunk = 128;
  int N_basic = 16;

  while((c = getopt_long(argc, argv, short_options, long_options, NULL) != -1)) {
    switch(c) {
    case 'c':
      N_chunk = strtol(optarg, NULL, 10);
      break;

    case 'b':
      N_basic = strtol(optarg, NULL, 10);
      break;

    case 'h':
      printf("Usage:\n");
      printf("\n");
      printf("{ --N_chunk | -c } N   The chunk size (currently %d)\n", N_chunk);
      printf("{ --N_basic | -b } N   The basic submatrix size (currently %d)\n", N_basic);
      printf("{ --help | -h }        This help.\n");
      exit(0);
      break;

    default:
      fprintf(stderr, "unknown option\n");
      exit(-1);
    }
  }
  return 0;
}
