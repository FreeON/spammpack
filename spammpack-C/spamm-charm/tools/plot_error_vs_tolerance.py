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
  end_of_run = True
  for line in fd:
    result = re.compile("read BCSR matrix.+[/]([^/]+).F_DIIS").search(line)
    if result:
      if not end_of_run:
        raise Exception("FIXME")
      end_of_run = False
      filename = result.group(1)
      logging.debug("new run, resetting fields, read " + filename)
      logging.debug("line = {:s}".format(line.rstrip()))
      if not filename in dataset:
        dataset[filename] = {}
      tolerance = None
      total_energy = None
      total_energy_D = None
      total_time = None
      multiply_time = 0
      complexity = None
      timers = {}

    result = re.compile("tolerance = ([0-9.e+-]+)$").search(line)
    if result:
      tolerance = float(result.group(1))
      logging.debug("read tolerance = {:e}".format(tolerance))

    result = re.compile("total energy, trace\(F\*P\)\s+=\s+(.*)$").search(line)
    if result:
      total_energy = float(result.group(1))
      logging.debug("read total energy = {:e}".format(total_energy))

    result = re.compile("total energy, trace\(F\*D\)\s+=\s+(.*)$").search(line)
    if result:
      total_energy_D = float(result.group(1))
      logging.debug("read total energy (D) = {:e}".format(total_energy_D))

    result = re.compile("total time.*:\s+(.*) seconds$").search(line)
    if result:
      total_time = float(result.group(1))
      logging.debug("read total time = {:e}".format(total_time))

    result = re.compile("complexity\s+=\s+(.*)$").search(line)
    if result:
      complexity = float(result.group(1))
      logging.debug("read complexity = {:e}".format(complexity))

    result = re.compile(
        "iteration\s+([0-9]+), multiply:\s([0-9.e+-]+) seconds"
        ).search(line)
    if result:
      iteration = int(result.group(1))
      if iteration in timers:
        raise Exception("FIXME")
      timers["multiply"] = float(result.group(2))
      multiply_time += float(result.group(2))
      logging.debug("read iteration {:d}".format(iteration))

    result = re.compile("end of spamm-charm").search(line)
    if result:
      logging.debug("end of run")
      end_of_run = True

      if tolerance in dataset[filename]:
        raise Exception("FIXME")

      dataset[filename][tolerance] = {}
      dataset[filename][tolerance]["timers"] = timers
      dataset[filename][tolerance]["energy"] = total_energy
      dataset[filename][tolerance]["energy (D)"] = total_energy_D
      dataset[filename][tolerance]["time"] = total_time
      dataset[filename][tolerance]["complexity"] = complexity
      dataset[filename][tolerance]["t_multiply"] = multiply_time

      if tolerance == 0 and total_energy_D:
        E_0 = dataset[filename][0]["energy"]
        BCSR_E = dataset[filename][0]["energy (D)"]
        dataset[filename][0]["BCSR energy error"] = abs((BCSR_E-E_0)/E_0)

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

  print("{:>10} {:>12} {:>13} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12} {:>12}".format(
    "filename",
    "tolerance",
    "energy",
    "rel. error",
    "BCSR error",
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
          + "{:e} ".format(dataset[filename][0]["BCSR energy error"])
          + "{:e} ".format(dataset[filename][tolerance]["t_multiply"])
          + "{:e} ".format(dataset[filename][tolerance]["time"])
          + "{:e} ".format(dataset[filename][tolerance]["complexity"])
          + "{:e} ".format(dataset[filename][tolerance]["relative t_multiply"])
          + "{:e} ".format(dataset[filename][tolerance]["relative time"])
          + "{:e} ".format(dataset[filename][tolerance]["relative complexity"])
        )

if __name__ == "__main__":
  main()
