#!/usr/bin/env python

class scaling_data:

  def __init__ (self):
    self.data = []

  def append (self):
    self.data.append({})

  def set_lambda (self, l):
    self.data[-1]["lambda"] = l

  def set_complexity (self, c):
    self.data[-1]["complexity"] = c

  def set_threads (self, P):
    self.data[-1]["threads"] = P

  def set_tolerance (self, t):
    self.data[-1]["tolerance"] = t

  def set_walltime (self, T):
    self.data[-1]["walltime"] = T

  def __str__ (self):
    return str(self.data)

  def get_complexity (self):
    c = []
    for i in self.data:
      if not i["complexity"] in c:
        c.append(i["complexity"])
    return c

  def get_threads (self):
    t = []
    for i in self.data:
      if not i["threads"] in t:
        t.append(i["threads"])
    return t

  def get_walltime (complexity = None, threads = None):
    print(complexity)
    print(threads)
    result = []
    for i in self.data:
      next_result = i
      if complexity and i["complexity"] != complexity:
        next_result = None
      if threads and i["threads"] != threads:
        next_result = None
      if next_result:
        result.append(next_result)
    return result

def main ():
  import argparse
  import re

  parser = argparse.ArgumentParser()

  parser.add_argument("FILE",
      help = "The output file of the scaling script")

  options = parser.parse_args()

  data = scaling_data()
  fd = open(options.FILE)
  for line in fd:
    if re.compile("running:").search(line):
      data.append()

    result = re.compile("lambda = (.*)$").search(line)
    if result:
      data.set_lambda(float(result.group(1)))

    result = re.compile("complexity ratio = (.*)$").search(line)
    if result:
      data.set_complexity(float(result.group(1)))

    result = re.compile("using ([0-9]*) OpenMP.*tolerance (.*), (.*) seconds").search(line)
    if result:
      data.set_threads(int(result.group(1)))
      data.set_tolerance(float(result.group(2)))
      data.set_walltime(float(result.group(3)))

  complexity_values = data.get_complexity()
  thread_values = data.get_threads()

  print(complexity_values)

  for t in thread_values:
    for c in complexity_values:
      print(c)
      walltime = data.get_walltime(complexity = c, threads = t)
      print("{:d} {:f}".format(c, walltime))

if __name__ == "__main__":
  main()
