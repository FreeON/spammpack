#!/usr/bin/env python

def main ():
  import argparse
  import logging
  import re

  parser = argparse.ArgumentParser(
      description = "This script parses the results of scaling tests "
      + "and collects relevant output data into one table for further "
      + "processing."
  )

  parser.add_argument(
      "FILE",
      help = "The scaling tests output file(s)",
      nargs = "+"
      )

  parser.add_argument(
      "--log",
      help = "Set the log level",
      choices = [ "INFO", "DEBUG" ],
      default = "INFO"
      )

  options = parser.parse_args()

  # Flatten filename list.
  while True:
    temp = []
    is_flat = True
    while len(options.FILE) > 0:
      i = options.FILE.pop(0)
      if type(i) == type([]):
        is_flat = False
        for j in i:
          temp.append(j)
      else:
        temp.append(i)
    for i in temp:
      options.FILE.append(i)
    if is_flat:
      break

  # Set the log level.
  logging.basicConfig(level = getattr(logging, options.log.upper(), None))

  logging.debug("starting...")
  logging.debug("parsing {:s}".format(options.FILE))

  dataset = {}

  for filename in options.FILE:
    fd = open(filename)
    end_of_run = True
    for line in fd:
      result = re.compile("read BCSR matrix.+[/]([^/]+).F_DIIS").search(line)
      if result:
        if not end_of_run:
          raise Exception("FIXME")
        end_of_run = False
        short_filename = result.group(1)
        logging.debug("new run, resetting fields, read " + short_filename)
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
          dataset[N]["filename"] = short_filename
          dataset[N]["data"] = {}

        if tolerance in dataset[N]["data"]:
          logging.warn("tolerance {:e} already in dataset".format(tolerance))
          logging.warn("while reading " + filename)
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
          logging.debug(
              "tolerance = {:e}, ".format(tolerance)
              + "E = {:e}, ".format(BCSR_E)
              + "E_0 = {:e}, ".format(E_0)
              + "BCSR error = {:e}".format(dataset[N]["data"][0]["BCSR energy error"])
              )

        if 0 in dataset[N]["data"]:
          E_0 = dataset[N]["data"][0]["energy"]
          E = dataset[N]["data"][tolerance]["energy"]
          dataset[N]["data"][tolerance]["energy error"] = abs((E-E_0)/E_0)
          logging.debug(
              "tolerance = {:e}, ".format(tolerance)
              + "E = {:e}, ".format(E)
              + "E_0 = {:e}, ".format(E_0)
              + "error = {:e}".format(dataset[N]["data"][tolerance]["energy error"])
              )

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

  Ns = sorted(dataset.keys())
  for N in Ns:
    tolerances = sorted(dataset[N]["data"].keys())

    BCSR_energy_error = 0
    if 0 in tolerances:
      if "BCSR energy error" in dataset[N]["data"][0]:
        BCSR_energy_error = dataset[N]["data"][0]["BCSR energy error"]

    print("# " + dataset[N]["filename"])
    print("# {:>3} {:>10} {:>10} {:>11} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}".format(
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

    for tolerance in tolerances:
      energy_error = 0
      if "energy error" in dataset[N]["data"][tolerance]:
        energy_error = dataset[N]["data"][tolerance]["energy error"]

      relative_complexity = 0
      if "relative complexity" in dataset[N]["data"][tolerance]:
        relative_complexity = dataset[N]["data"][tolerance]["relative complexity"]

      relative_t_multiply = 0
      if "relative t_multiply" in dataset[N]["data"][tolerance]:
        relative_t_multiply = dataset[N]["data"][tolerance]["relative t_multiply"]

      relative_time = 0
      if "relative time" in dataset[N]["data"][tolerance]:
        relative_time = dataset[N]["data"][tolerance]["relative time"]

      print(
          "{:5d} ".format(N)
          + "{:>10} ".format(dataset[N]["filename"])
          + "{:1.4e} ".format(tolerance)
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["energy"])
          + "{:1.4e} ".format(energy_error)
          + "{:1.4e} ".format(BCSR_energy_error)
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["t_multiply"])
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["t_symbolic_multiply"])
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["time"])
          + "{:1.4e} ".format(dataset[N]["data"][tolerance]["complexity"])
          + "{:1.4e} ".format(relative_t_multiply)
          + "{:1.4e} ".format(relative_time)
          + "{:1.4e}".format(relative_complexity)
        )

    print("")
    print("")

if __name__ == "__main__":
  main()
