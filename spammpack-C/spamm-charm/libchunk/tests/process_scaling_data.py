#!/usr/bin/env python

class scaling_data:

  def __init__ (self, filename):
    import re

    self.data = []

    fd = open(filename)
    for line in fd:
      if re.compile("running:").search(line):
        self.append()

      result = re.compile("lambda = (.*)$").search(line)
      if result:
        self.set_lambda(float(result.group(1)))

      result = re.compile("complexity ratio = (.*)$").search(line)
      if result:
        self.set_complexity(float(result.group(1)))

      result = re.compile("using ([0-9]*) OpenMP.*tolerance (.*), (.*) seconds").search(line)
      if result:
        self.set_threads(int(result.group(1)))
        self.set_tolerance(float(result.group(2)))
        self.set_walltime(float(result.group(3)))

    fd.close()

  def info (self):
    info  = "complexity: " + str(self.get_complexity()) + "\n"
    info += "thread:     " + str(self.get_threads())
    return info

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
    result = ""
    for i in self.data:
      result += str(i) + "\n"
    return result

  def get_complexity (self):
    c = []
    for i in self.data:
      if not i["complexity"] in c:
        c.append(i["complexity"])
    return sorted(c, reverse = True)

  def get_threads (self):
    t = []
    for i in self.data:
      if not i["threads"] in t:
        t.append(i["threads"])
    return sorted(t)

  def get_walltime (self, complexity = None, threads = None):
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

def flatten_list (l):
  while True:
    temp = []
    isFlat = True
    for i in l:
      if type(i) == type([]):
        isFlat = False
        for j in i:
          temp.append(j)
      else:
        temp.append(i)
    l = temp
    if isFlat:
      break
  return l

def main ():
  import argparse
  import matplotlib.pyplot as plt

  parser = argparse.ArgumentParser()

  parser.add_argument("FILE",
      help = "The output file of the scaling script")

  parser.add_argument("--thread",
      metavar = "N",
      help = "Plot only results for N threads",
      nargs = "+",
      action = "append")

  parser.add_argument("--complexity",
      metavar = "C",
      help = "Plot only results for complexity C",
      nargs = "+",
      action = "append")

  parser.add_argument("--type",
      help = "The plot type (default %(default)s)",
      choices = [
        "speedup_complexity",
        "speedup_thread"
        ],
      default = "speedup_complexity")

  options = parser.parse_args()

  if options.thread:
    options.thread = flatten_list(options.thread)
    print("plotting only threads " + str(options.thread))

  if options.complexity:
    options.complexity = flatten_list(options.complexity)
    print("plotting only complexity " + str(options.complexity))

  data = scaling_data(options.FILE)
  print(data.info())

  if options.complexity:
    complexity_values = sorted(
        [ float(i) for i in options.complexity ],
        reverse = True
        )
  else:
    complexity_values = data.get_complexity()

  if options.thread:
    thread_values = sorted([ int(i) for i in options.thread ])
  else:
    thread_values = data.get_threads()

  # Plot walltime vs. complexity.
  plt.figure()

  if options.type == "speedup_complexity":
    for t in thread_values:
      walltime = []
      for c in complexity_values:
        query = data.get_walltime(complexity = c, threads = t)
        if len(query) != 1:
          raise Exception("can not find result for {:d} threads".format(t))
        walltime.append(query[0]["walltime"])
      plt.plot(
          complexity_values,
          [ walltime[0]/i for i in walltime ],
          linestyle = "-",
          marker = "o",
          label = "{:d} threads".format(t)
          )

    plt.plot(
        complexity_values,
        [ 1/i for i in complexity_values ],
        label = "ideal"
        )

    plt.gca().invert_xaxis()
    plt.legend()
    plt.xlabel("complexity")
    plt.ylabel("speedup vs. dense")
    plt.show()

  if options.type == "speedup_thread":
    for c in complexity_values:
      walltime = []
      for t in thread_values:
        query = data.get_walltime(complexity = c, threads = t)
        if len(query) != 1:
          raise Exception("can not find result for {:d} threads".format(t))
        walltime.append(query[0]["walltime"])
      plt.plot(
          thread_values,
          [ walltime[0]/i for i in walltime ],
          linestyle = "-",
          marker = "o",
          label = "complexity {:1.3f}".format(c)
          )

    plt.plot(
        thread_values,
        thread_values,
        label = "ideal"
        )

    plt.legend()
    plt.xlabel("threads")
    plt.ylabel("speedup vs. dense")
    plt.show()

  else:
    raise Exception("unknown plot type")

if __name__ == "__main__":
  main()
