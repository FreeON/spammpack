#!/usr/bin/env python

import random

for i in range(10):
  for j in range(10):
    print "%d %d % 1.14e" % (i+1, j+1, random.random())
