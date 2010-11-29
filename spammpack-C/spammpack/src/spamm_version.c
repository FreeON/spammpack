#include "config.h"
#include "spamm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *
spamm_version ()
{
  char *string = (char *) malloc(sizeof(char)*strlen(VERSION)+1+strlen(COMMIT_TAG)+1);

  sprintf(string, "%s:%s", VERSION, COMMIT_TAG);
  return string;
}
