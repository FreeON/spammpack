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
      N = None
      tolerance = None
      total_energy = None
      total_energy_D = None
      total_time = None
      multiply_time = 0
      symbolic_multiply_time = 0
      complexity = None
      timers = {}

    result = re.compile("name = F, N = ([0-9]+)").search(line)
    if result:
      N = int(result.group(1))
      logging.debug("read N = {:d}".format(N))

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
        "iteration\s+([0-9]+), "
        + "multiply:\s+([0-9.e+-]+) seconds, "
        + "symbolic multiply:\s+([0-9.e+-]+) seconds, "
        ).search(line)
    if result:
      iteration = int(result.group(1))
      if iteration in timers:
        raise Exception("FIXME")
      timers["multiply"] = float(result.group(2))
      multiply_time += float(result.group(2))
      symbolic_multiply_time += float(result.group(3))
      logging.debug("read iteration {:d}".format(iteration))

    result = re.compile("end of spamm-charm").search(line)
    if result:
      logging.debug("end of run")
      end_of_run = True

      if not N in dataset:
        dataset[N] = {}
        dataset[N]["filename"] = filename
        dataset[N]["data"] = {}

      if tolerance in dataset[N]["data"]:
        raise Exception("FIXME")

      dataset[N]["data"][tolerance] = {}
      dataset[N]["data"][tolerance]["timers"] = timers
      dataset[N]["data"][tolerance]["energy"] = total_energy
      dataset[N]["data"][tolerance]["energy (D)"] = total_energy_D
      dataset[N]["data"][tolerance]["time"] = total_time
      dataset[N]["data"][tolerance]["complexity"] = complexity
      dataset[N]["data"][tolerance]["t_multiply"] = multiply_time
      dataset[N]["data"][tolerance]["t_symbolic_multiply"] = symbolic_multiply_time

      if tolerance == 0 and total_energy_D:
        E_0 = dataset[N]["data"][0]["energy"]
        BCSR_E = dataset[N]["data"][0]["energy (D)"]
        dataset[N]["data"][0]["BCSR energy error"] = abs((BCSR_E-E_0)/E_0)

      if 0 in dataset[N]["data"]:
        E_0 = dataset[N]["data"][0]["energy"]
        E = dataset[N]["data"][tolerance]["energy"]
        dataset[N]["data"][tolerance]["energy error"] = abs((E-E_0)/E_0)

        t_0 = dataset[N]["data"][0]["time"]
        t = dataset[N]["data"][tolerance]["time"]
        dataset[N]["data"][tolerance]["relative time"] = t/t_0

        t_0 = dataset[N]["data"][0]["t_multiply"]
        t = dataset[N]["data"][tolerance]["t_multiply"]
        dataset[N]["data"][tolerance]["relative t_multiply"] = t/t_0

        c_0 = dataset[N]["data"][0]["complexity"]
        c = dataset[N]["data"][tolerance]["complexity"]
        dataset[N]["data"][tolerance]["relative complexity"] = c/c_0

  fd.close()

  print("{:>5} {:>10} {:>10} {:>11} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}".format(
    "N",
    "filename",
    "tolerance",
    "energy",
    "SpAMM_err",
    "BCSR_err",
    "t_mul",
    "t_sym",
    "t_tot",
    "compl.",
    "d(t_mul)",
    "d(t_tot)",
    "d(compl.)"
    ))

  Ns = sorted(dataset.keys())
  for N in Ns:
    tolerances = sorted(dataset[N]["data"].keys())
    for tolerance in tolerances:
      print(
          "{:5d} ".format(N)
          + "{:>10} ".format(dataset[N]["filename"])
          + "{:1.4e} ".format(tolerance)
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["energy"])
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["energy error"])
          + "{:1.4e} ".format(dataset[N]["data"][0]["BCSR energy error"])
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["t_multiply"])
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["t_symbolic_multiply"])
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["time"])
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["complexity"])
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["relative t_multiply"])
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["relative time"])
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["relative complexity"])
        )

if __name__ == "__main__":
  main()
