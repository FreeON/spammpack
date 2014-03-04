#!/usr/bin/env python

def main ():
  import argparse
  import logging
  import re

  parser = argparse.ArgumentParser()

  parser.add_argument(
      "FILE",
      help = "The scaling output file"
      )

  parser.add_argument(
      "--log",
      help = "Set the log level",
      choices = [ "INFO", "DEBUG" ],
      default = "INFO"
      )

  options = parser.parse_args()

  # Set the log level.
  logging.basicConfig(level = getattr(logging, options.log.upper(), None))

  logging.debug("starting...")

  fd = open(options.FILE)
  dataset = {}
  for line in fd:
    result = re.compile("read BCSR matrix.+[/]([^/]+).F_DIIS").search(line)
    if result:
      logging.debug("new run, resetting fields")
      filename = result.group(1)
      if not filename in dataset:
        dataset[filename] = {}
      tolerance = None
      total_energy = None
      total_time = None
      multiply_time = 0
      complexity = None
      timers = {}

    result = re.compile("tolerance = ([0-9.e+-]+)$").search(line)
    if result:
      tolerance = float(result.group(1))
      logging.debug("read tolerance = {:e}".format(tolerance))

    result = re.compile("total energy.*=\s+(.*)$").search(line)
    if result:
      total_energy = float(result.group(1))

    result = re.compile("total time.*:\s+(.*) seconds$").search(line)
    if result:
      total_time = float(result.group(1))

    result = re.compile("complexity\s+=\s+(.*)$").search(line)
    if result:
      complexity = float(result.group(1))

    result = re.compile(
        "iteration\s([0-9]+), multiply:\s([0-9.e+-]+) seconds"
        ).search(line)
    if result:
      iteration = int(result.group(1))
      if iteration in timers:
        raise Exception("FIXME")
      timers["multiply"] = float(result.group(2))
      multiply_time += float(result.group(2))

    result = re.compile("Program finished").search(line)
    if result:
      if tolerance in dataset[filename]:
        raise Exception("FIXME")
      dataset[filename][tolerance] = {}
      dataset[filename][tolerance]["timers"] = timers
      dataset[filename][tolerance]["energy"] = total_energy
      dataset[filename][tolerance]["time"] = total_time
      dataset[filename][tolerance]["complexity"] = complexity
      dataset[filename][tolerance]["t_multiply"] = multiply_time

      if 0 in dataset[filename]:
        E_0 = dataset[filename][0]["energy"]
        E = dataset[filename][tolerance]["energy"]
        dataset[filename][tolerance]["energy error"] = abs((E-E_0)/E_0)

        t_0 = dataset[filename][0]["time"]
        t = dataset[filename][tolerance]["time"]
        dataset[filename][tolerance]["relative time"] = t/t_0

        t_0 = dataset[filename][0]["t_multiply"]
        t = dataset[filename][tolerance]["t_multiply"]
        dataset[filename][tolerance]["relative t_multiply"] = t/t_0

        c_0 = dataset[filename][0]["complexity"]
        c = dataset[filename][tolerance]["complexity"]
        dataset[filename][tolerance]["relative complexity"] = c/c_0

  fd.close()

  print("{:>10} {:>12} {:>13} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12}".format(
    "filename",
    "tolerance",
    "energy",
    "rel. error",
    "t_multiply",
    "total time",
    "complexity",
    "rel. t_mult.",
    "rel. time",
    "rel. compl."
    ))

  filenames = sorted(dataset.keys())
  for filename in filenames:
    tolerances = sorted(dataset[filename].keys())
    for tolerance in tolerances:
      print(
          "{:>10} ".format(filename)
          + "{:e} ".format(tolerance)
          + "{:e} ".format(dataset[filename][tolerance]["energy"])
          + "{:e} ".format(dataset[filename][tolerance]["energy error"])
          + "{:e} ".format(dataset[filename][tolerance]["t_multiply"])
          + "{:e} ".format(dataset[filename][tolerance]["time"])
          + "{:e} ".format(dataset[filename][tolerance]["complexity"])
          + "{:e} ".format(dataset[filename][tolerance]["relative t_multiply"])
          + "{:e} ".format(dataset[filename][tolerance]["relative time"])
          + "{:e} ".format(dataset[filename][tolerance]["relative complexity"])
        )

if __name__ == "__main__":
  main()
