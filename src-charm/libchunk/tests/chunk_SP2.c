#include <getopt.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#define ABORT(format, ...) print_log("ERROR", __FILE__, __LINE__, __func__, format, ##__VA_ARGS__); exit(-1)
#define INFO(format, ...) print_log("INFO", __FILE__, __LINE__, __func__, format, ##__VA_ARGS__)

#define FORMAT_LENGTH 1000

#include "bcsr.h"

void
print_log (const char *const tag,
    const char *const filename,
    const int line_number,
    const char *const function_name,
    const char *const format,
    ...)
{
  va_list ap;
  char new_format[FORMAT_LENGTH];

  snprintf(new_format, FORMAT_LENGTH, "[%s - %s:%d (%s)] %s",
      tag, filename, line_number, function_name, format);

  va_start(ap, format);
  vprintf(new_format, ap);
  va_end(ap);
}

void
print_help_and_exit (void)
{
  printf("Usage: chunk_SP2 [options]\n");
  printf("\n");
  printf("{ -h | --help }     This help\n");
  exit(0);
}

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
      case 'h':
        print_help_and_exit();
        break;

      default:
        printf("illegal command line argument\n");
        exit(-1);
        break;
    };
  }

  //if(F_filename == NULL)
  //{
  //  ABORT("missing Fockian matrix file\n");
  //}

  return 0;
}
