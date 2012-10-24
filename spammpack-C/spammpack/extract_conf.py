#!/usr/bin/env python

import argparse
import re

parser = argparse.ArgumentParser()

parser.add_argument("CONFFILE",
    help = "The config file (default = %(default)s)",
    default = "config.h")

options = parser.parse_args()

fd = open(options.CONFFILE)
lines = fd.readlines()
fd.close()

print("!> @file")
print()
print("#ifndef __SPAMM_CONFIG_H")
print("#define __SPAMM_CONFIG_H")
print()
print("! This file contains the configuration options in a Fortran")
print("! friendly format.")
print()

for i in range(len(lines)):
  result = re.compile("^#define").search(lines[i])
  if result:
    # Look for preceeding comment.
    j = i-1
    while j >= 0:
      result = re.compile("^/\*").search(lines[j])
      if result:
        # Doxygenify the comment.
        while j <= i:
          doxyfied_line = lines[j].rstrip()
          doxyfied_line = re.sub("^/\* +", "!> ", doxyfied_line)
          doxyfied_line = re.sub("^ +", "!> ", doxyfied_line)
          doxyfied_line = re.sub("\*/.*$", "", doxyfied_line)
          if doxyfied_line != "!> ":
            print(doxyfied_line)
          j += 1
        break
      j -= 1
    print()
  i += 1

print("#endif")
