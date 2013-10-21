#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy
import os
import re
import signal
import sys

def make_graph ():
  for PPV in timing:
    np = []
    for i in timing[PPV]:
      np.append(i)
    np.sort()

    t = []
    np_dead = []
    if iteration < 0:
      additional_label = "all"
      for i in np:
        t.append(timing[PPV][i]["t_total"])
    else:
      additional_label = "iteration {:d}".format(iteration)
      for i in np:
        if iteration <= len(timing[PPV][i]["t_multiply"]):
          t.append(
              timing[PPV][i]["t_multiply"][iteration-1]
              +timing[PPV][i]["t_add"][iteration-1]
              +timing[PPV][i]["t_setEq"][iteration-1])
        else:
          np_dead.append(i)

    for i in np_dead:
      np.remove(i)

    plt.loglog(np, t, '-o', label = "PPV{:03d} ({:s})".format(
      PPV, additional_label), hold = True)
    plt.loglog(
        [ np[0], np[-1] ],
        [ t[0], t[0]/np[-1] ],
        '-',
        hold = True)

parser = argparse.ArgumentParser()

parser.add_argument("OUTFILE",
    help = "The output file of a spamm SP2 calculation",
    nargs = "+")

parser.add_argument("--iteration",
    metavar = "N",
    help = "Plot iteration N instead of total time",
    type = int,
    nargs = "+",
    action = "append")

parser.add_argument("--one-graph",
    help = "Plot all iterations on one graph",
    action = "store_true",
    default = False)

parser.add_argument("--debug",
    help = "Print debug information",
    action = "store_true",
    default = False)

options = parser.parse_args()

while True:
  temp = []
  flattened = True
  for i in options.iteration:
    if type(i) != type(1):
      flattened = False
      for j in i:
        temp.append(j)
    else:
      temp.append(i)
  options.iteration = temp
  if flattened:
    break

timing = {}

for f in options.OUTFILE:
  fd = open(f)
  np = 0
  PPV = 0
  t_multiply = []
  t_add = []
  t_setEq = []
  for line in fd:
    if np == 0:
      result = re.compile("Running on\s+([0-9]+)\s+unique").search(line)
      if result:
        np = int(result.group(1))

    if PPV == 0:
      result = re.compile("PP([0-9]+).OrthoF").search(line)
      if result:
        PPV = int(result.group(1))
        if not PPV in timing:
          timing[PPV] = {}

    result = re.compile(
        "iteration\s+([0-9]+), "
        + "multiply: ([0-9.]+) seconds, "
        + "setEq: ([0-9.]+) seconds, "
        + "add: ([0-9.]+) seconds").search(line)
    if result:
      iteration = int(result.group(1))
      t_multiply.append(float(result.group(2)))
      t_add.append(float(result.group(3)))
      t_setEq.append(float(result.group(4)))

      if iteration != len(t_multiply):
        raise Exception("strange...")

  t_multiply = numpy.array(t_multiply)
  t_add = numpy.array(t_add)
  t_setEq = numpy.array(t_setEq)

  timing[PPV][np] = {
      "t_total" : numpy.sum(t_multiply)+numpy.sum(t_add)+numpy.sum(t_setEq),
      "t_multiply" : t_multiply,
      "t_add" : t_add,
      "t_setEq" : t_setEq
      }

  if options.debug:
    print("file {:s}, PPV{:03d}, np = {:d}, {:d} iterations".format(
      f, PPV, np, len(timing[PPV][np]["t_multiply"])))

pid = []
for iteration in options.iteration:
  if not options.one_graph:
    pid.append(os.fork())

    if(pid[-1]):
      pass

    else:
      plt.figure(iteration)
      make_graph()
      plt.legend()
      plt.show()

  else:
    make_graph()


if options.one_graph:
  plt.legend()
  plt.show()

else:
  print("done, please press Enter to close all windows")
  sys.stdin.readline()

  for i in pid:
    os.kill(i, signal.SIGHUP)
