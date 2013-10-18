#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy
import re

parser = argparse.ArgumentParser()

parser.add_argument("OUTFILE",
    help = "The output file of a spamm SP2 calculation",
    nargs = "+")

options = parser.parse_args()

timing = {}

for f in options.OUTFILE:
  fd = open(f)
  np = 0
  PPV = 0
  multiply = []
  add = []
  setEq = []
  for line in fd:
    if np == 0:
      result = re.compile("Reported:.*\(out of ([0-9]+)\)").search(line)
      if result:
        np = int(result.group(1))
    if PPV == 0:
      result = re.compile("PP([0-9]+).OrthoF").search(line)
      if result:
        PPV = int(result.group(1))
    result = re.compile("iteration\s+([0-9]+), "
        + "multiply: ([0-9.]+) seconds, "
        + "setEq: ([0-9.]+) seconds, "
        + "add: ([0-9.]+) seconds").search(line)
    if result:
      if len(multiply) != int(result.group(1))-1:
        raise Exception("illegal iteration, len(multiply) = {:d}, iteration = {:d}".format(
          len(multiply), int(result.group(1))))
      multiply.append(float(result.group(2)))
      add.append(float(result.group(3)))
      setEq.append(float(result.group(4)))

  t_multiply = numpy.array(multiply)
  t_add = numpy.array(add)
  t_setEq = numpy.array(setEq)

  if not PPV in timing:
    timing[PPV] = { "np" : [], "t_total" : [] }

  timing[PPV]["np"].append(np)
  timing[PPV]["t_total"].append(numpy.sum(t_multiply)
      + numpy.sum(t_add) + numpy.sum(t_setEq))

for PPV in timing:
  print("adding PPV{:03d}".format(PPV))
  plt.loglog(timing[PPV]["np"],
      timing[PPV]["t_total"],
      '-o',
      label = "PPV{:03d}".format(PPV))
  plt.hold()
  plt.loglog([ timing[PPV]["np"][0], timing[PPV]["np"][-1] ],
      [ timing[PPV]["t_total"][0], timing[PPV]["t_total"][0]/timing[PPV]["np"][-1] ],
      '-',
      label = "PPV{:03d} ideal".format(PPV))
  plt.hold()

plt.legend()
plt.show()
