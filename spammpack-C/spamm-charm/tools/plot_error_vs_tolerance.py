#!/usr/bin/env python

def main ():
  import argparse
  import re

  parser = argparse.ArgumentParser()

  parser.add_argument(
      "FILE",
      help = "The scaling output file"
      )

  options = parser.parse_args()

  fd = open(options.FILE)
  dataset = {}
  for line in fd:
    result = re.compile("read BCSR matrix.+[/]([^/]+).F_DIIS").search(line)
    if result:
      filename = result.group(1)
      if not filename in dataset:
        dataset[filename] = {}
      tolerance = None
      total_energy = None
      total_time = None
      complexity = None

    result = re.compile("tolerance\s*=\s+(.*)$").search(line)
    if result:
      tolerance = float(result.group(1))

    result = re.compile("total energy.*=\s+(.*)$").search(line)
    if result:
      total_energy = float(result.group(1))

    result = re.compile("total time.*:\s+(.*) seconds$").search(line)
    if result:
      total_time = float(result.group(1))

    result = re.compile("complexity\s+=\s+(.*)$").search(line)
    if result:
      complexity = int(result.group(1))

    result = re.compile("Program finished").search(line)
    if result:
      if tolerance in dataset[filename]:
        raise Exception("FIXME")
      dataset[filename][tolerance] = {}
      dataset[filename][tolerance]["energy"] = total_energy
      dataset[filename][tolerance]["time"] = total_time
      dataset[filename][tolerance]["complexity"] = complexity

      if 0 in dataset[filename]:
        E_0 = dataset[filename][0]["energy"]
        E = dataset[filename][tolerance]["energy"]
        dataset[filename][tolerance]["energy error"] = abs((E-E_0)/E_0)

        t_0 = dataset[filename][0]["time"]
        t = dataset[filename][tolerance]["time"]
        dataset[filename][tolerance]["relative time"] = t/t_0

        c_0 = dataset[filename][0]["complexity"]
        c = dataset[filename][tolerance]["complexity"]
        dataset[filename][tolerance]["relative complexity"] = c/c_0

  fd.close()

  print("{:>10} {:>12} {:>13} {:>12} {:>12} {:>12} {:>12} {:>12}".format(
    "filename",
    "tolerance",
    "energy",
    "rel. error",
    "timer",
    "rel. time",
    "complexity",
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
          + "{:e} ".format(dataset[filename][tolerance]["time"])
          + "{:e} ".format(dataset[filename][tolerance]["relative time"])
          + "{:12d} ".format(dataset[filename][tolerance]["complexity"])
          + "{:e}".format(dataset[filename][tolerance]["relative complexity"])
        )

if __name__ == "__main__":
  main()
