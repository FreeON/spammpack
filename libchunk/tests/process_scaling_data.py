#!/usr/bin/env python

class scaling_data:

  def __init__ (self, filename):
    import re

    self.data = []
    self.N_chunk = 0
    self.N_basic = 0

    fd = open(filename)
    for line in fd:
      if self.N_chunk == 0:
        result = re.compile("^running.* -N ([0-9]+) -b ([0-9]*) ").search(line)
        if result:
          self.N_chunk = int(result.group(1))
          self.N_basic = int(result.group(2))

      if re.compile("running:").search(line):
        self.append()

      result = re.compile("lambda = (.*)$").search(line)
      if result:
        self.set_lambda(float(result.group(1)))

      result = re.compile("complexity ratio = (.*)$").search(line)
      if result:
        self.set_complexity(float(result.group(1)))

      result = re.compile("done multiplying using ([0-9]*) OpenMP.*tolerance (.*), (.*) seconds").search(line)
      if result:
        self.set_threads(int(result.group(1)))
        self.set_tolerance(float(result.group(2)))
        self.set_dense(False)
        self.set_walltime(float(result.group(3)))

      result = re.compile("done multiplying dense using ([0-9]*) OpenMP.* (.*) seconds").search(line)
      if result:
        self.set_threads(int(result.group(1)))
        self.set_dense(True)
        self.set_tolerance(0.0)
        self.set_walltime(float(result.group(2)))

    fd.close()

  def info (self):
    info = ""
    info += "N = {:d}, N_basic = {:d}\n".format(self.N_chunk, self.N_basic)
    info += "complexity: " + str(self.get_complexity()) + "\n"
    info += "thread:     " + str(self.get_threads())
    return info

  def append (self):
    self.data.append({})

  def set_dense (self, isDense):
    self.data[-1]["isDense"] = isDense

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

  def get_tolerance (self):
    tau = []
    for i in self.data:
      if not i["tolerance"] in tau:
        tau.append(i["tolerance"])
    return sorted(tau)

  def get_threads (self):
    t = []
    for i in self.data:
      if not i["threads"] in t:
        t.append(i["threads"])
    return sorted(t)

  def get_record (
      self,
      isDense = False,
      complexity = None,
      threads = None,
      tolerance = None):
    result = []
    for i in self.data:
      next_result = i
      if complexity and i["complexity"] != complexity:
        next_result = None
      if tolerance and i["tolerance"] != tolerance:
        next_result = None
      if threads and i["threads"] != threads:
        next_result = None
      if isDense != i["isDense"]:
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


def plot_walltime_vs_complexity (data, options):
  import matplotlib.pyplot as plt

  plt.figure(
      figsize = (
        options.width/100.0,
        options.height/100.0
        ),
      dpi = 100
      )

  complexity_values = data.get_complexity()
  if options.thread:
    thread_values = sorted([ int(i) for i in options.thread ])
  else:
    thread_values = data.get_threads()

  max_speedup = 1
  for t in thread_values:
    walltime = []
    for c in complexity_values:
      query = data.get_record(complexity = c, threads = t)
      if len(query) == 0:
        raise Exception("can not find result for "
        + "complexity {:1.3f} and {:d} threads".format(c, t))
      walltime.append(query[0]["walltime"])
    if max_speedup < max([ walltime[0]/i for i in walltime ]):
      max_speedup = max([ walltime[0]/i for i in walltime ])
    plt.loglog(
        complexity_values,
        [ walltime[0]/i for i in walltime ],
        linestyle = "-",
        marker = "o",
        label = "{:d} threads".format(t)
        )

  plt.loglog(
      complexity_values,
      [ 1/i for i in complexity_values ],
      color = "black",
      label = "ideal"
      )

  plt.grid(True)
  plt.xlim([min(complexity_values), max(complexity_values)])
  plt.ylim([1, max_speedup])
  plt.gca().invert_xaxis()
  plt.legend(loc = "upper left")
  plt.xlabel("complexity")
  plt.ylabel("parallel speedup")
  if not options.no_title:
    plt.title("N = {:d}, N_basic = {:d}".format(data.N_chunk, data.N_basic))

  if options.output:
    plt.savefig(options.output + "_walltime_vs_complexity.png")

def plot_walltime_vs_tolerance (data, options):
  import matplotlib.pyplot as plt

  plt.figure(
      figsize = (
        options.width/100.0,
        options.height/100.0
        ),
      dpi = 100
      )

  tolerance_values = data.get_tolerance()
  if options.thread:
    thread_values = sorted([ int(i) for i in options.thread ])
  else:
    thread_values = data.get_threads()

  max_speedup = 1
  for t in thread_values:
    walltime = []
    for c in tolerance_values:
      query = data.get_record(tolerance = c, threads = t)
      if len(query) == 0:
        raise Exception("can not find result for "
        + "tolerance {:e} and {:d} threads".format(c, t))
      walltime.append(query[0]["walltime"])
    if max_speedup < max([ walltime[0]/i for i in walltime ]):
      max_speedup = max([ walltime[0]/i for i in walltime ])
    plt.semilogx(
        tolerance_values,
        walltime,
        linestyle = "-",
        marker = "o",
        label = "{:d} threads".format(t)
        )

  plt.grid(True)
  plt.legend(loc = "upper right")
  plt.xlabel("tolerance")
  plt.ylabel("walltime [s]")
  if not options.no_title:
    plt.title("N = {:d}, N_basic = {:d}".format(data.N_chunk, data.N_basic))

  if options.output:
    plt.savefig(options.output + "_walltime_vs_tolerance.png")

def plot_walltime_vs_threads (data, options):
  import matplotlib.pyplot as plt

  plt.figure(
      figsize = (
        options.width/100.0,
        options.height/100.0
        ),
      dpi = 100
      )

  #if options.complexity:
  #  complexity_values = sorted(
  #      [ float(i) for i in options.complexity ],
  #      reverse = True
  #      )
  #else:
  #  complexity_values = data.get_complexity()
  tolerance_values = data.get_tolerance()
  thread_values = data.get_threads()

  #for c in complexity_values:
  for c in tolerance_values:
    if options.print:
      print("# tolerance = %e" % (c))
    walltime = []
    for t in thread_values:
      #query = data.get_record(complexity = c, threads = t)
      query = data.get_record(tolerance = c, threads = t)
      if len(query) == 0:
        raise Exception("can not find SpAMM result for {:d} threads".format(t))
      walltime.append(query[0]["walltime"])
    walltime_values = [ walltime[0]/i for i in walltime ]
    plt.plot(
        thread_values,
        walltime_values,
        linestyle = "-",
        marker = "o",
        label = "complexity {:1.3f}".format(c)
        )
    if options.print:
      for i in range(len(thread_values)):
        print("%3d %e" % (thread_values[i], walltime[i]))
      print("\n")

  walltime = []
  for t in thread_values:
    query = data.get_record(isDense = True, threads = t)
    if len(query) == 0:
      raise Exception("can not find dense result for {:d} threads".format(t))
    walltime.append(query[0]["walltime"])
  plt.plot(
      thread_values,
      [ walltime[0]/i for i in walltime ],
      linestyle = "-",
      marker = "*",
      label = "dense"
      )

  plt.plot(
      thread_values,
      thread_values,
      color = "black",
      label = "ideal"
      )

  plt.grid(True)
  plt.legend(loc = "upper left")
  plt.xlim([min(thread_values), max(thread_values)])
  plt.ylim([min(thread_values), max(thread_values)])
  plt.xlabel("threads")
  plt.ylabel("parallel speedup")
  if not options.no_title:
    plt.title("N = {:d}, N_basic = {:d}".format(data.N_chunk, data.N_basic))

  if options.output:
    plt.savefig(options.output + "_walltime_vs_threads.png")

def plot_walltime (data, options):
  import matplotlib.pyplot as plt

  figure, ax = plt.subplots()

  width = 0.3
  plt.bar(
      1-width/2,
      data.get_record(
        isDense = True,
        threads = 1
        )[0]["walltime"],
      width,
      color = "red"
      )
  plt.bar(
      2-width/2,
      data.get_record(
        complexity = 1,
        threads = 1
        )[0]["walltime"],
      width,
      color = "blue"
      )

  plt.bar(
      3-width/2,
      data.get_record(
        isDense = True,
        threads = max(data.get_threads())
        )[0]["walltime"],
      width,
      color = "red"
      )
  plt.bar(
      4-width/2,
      data.get_record(
        complexity = 1,
        threads = max(data.get_threads())
        )[0]["walltime"],
      width,
      color = "blue"
      )

  plt.xlabel('linear algebra library')
  plt.ylabel('walltime [s]')

  ax.set_xticks([ 1, 2, 3, 4 ])
  ax.set_xticklabels(
      (
        "MKL (serial)",
        "SpAMM (serial)",
        "MKL (%d threads)" % (max(data.get_threads())),
        "SpAMM (%d threads)" % (max(data.get_threads()))
        )
      )

  plt.xlim([ 0.5, 4.5 ])

  if options.output:
    plt.savefig(options.output + "_walltimg.png")

def plot_efficiency_vs_threads (data, options):
  import matplotlib.pyplot as plt

  plt.figure(
      figsize = (
        options.width/100.0,
        options.height/100.0
        ),
      dpi = 100
      )

  if options.complexity:
    complexity_values = sorted(
        [ float(i) for i in options.complexity ],
        reverse = True
        )
  else:
    complexity_values = data.get_complexity()
  thread_values = data.get_threads()

  for c in complexity_values:
    walltime = []
    for t in thread_values:
      query = data.get_record(complexity = c, threads = t)
      if len(query) == 0:
        raise Exception("can not find SpAMM result for {:d} threads".format(t))
      walltime.append(query[0]["walltime"])
    plt.plot(
        thread_values,
        [ 100*walltime[0]/walltime[i]/thread_values[i] for i in range(len(walltime)) ],
        linestyle = "-",
        marker = "o",
        label = "complexity {:1.3f}".format(c)
        )

  walltime = []
  for t in thread_values:
    query = data.get_record(isDense = True, threads = t)
    if len(query) == 0:
      raise Exception("can not find dense result for {:d} threads".format(t))
    walltime.append(query[0]["walltime"])
  plt.plot(
      thread_values,
      [ 100*walltime[0]/walltime[i]/thread_values[i] for i in range(len(walltime)) ],
      linestyle = "-",
      marker = "*",
      label = "dense"
      )

  plt.grid(True)
  plt.legend(loc = "lower left")
  plt.xlabel("threads")
  plt.ylabel("parallel efficiency")
  if not options.no_title:
    plt.title("N = {:d}, N_basic = {:d}".format(data.N_chunk, data.N_basic))

  if options.output:
    plt.savefig(options.output + "_efficiency_vs_threads.png")

def plot_efficiency_vs_complexity (data, options):
  import matplotlib.pyplot as plt

  plt.figure(
      figsize = (
        options.width/100.0,
        options.height/100.0
        ),
      dpi = 100
      )

  complexity_values = data.get_complexity()
  if options.thread:
    thread_values = sorted([ int(i) for i in options.thread ])
  else:
    thread_values = data.get_threads()

  if 1 in complexity_values:
    for t in thread_values:
      walltime = []
      walltime_1 = 0
      for i in range(len(complexity_values)):
        query = data.get_record(complexity = complexity_values[i], threads = t)
        if len(query) == 0:
          raise Exception("can not find result for "
          + "complexity {:1.3f} and {:d} threads".format(
            complexity_values[i], t
            )
          )
        walltime.append(query[0]["walltime"])
        if complexity_values[i] == 1:
          walltime_1 = walltime[-1]
      plt.plot(
          complexity_values,
          [ 100*complexity_values[i]*walltime_1/walltime[i] for i in range(len(walltime)) ],
          linestyle = "-",
          marker = "o",
          label = "{:d} threads".format(t)
          )

    plt.grid(True)
    plt.gca().invert_xaxis()
    plt.legend(loc = "lower left")
    plt.xlabel("complexity")
    plt.ylabel("complexity effciciency")
    if not options.no_title:
      plt.title("N = {:d}, N_basic = {:d}".format(data.N_chunk, data.N_basic))

    if options.output:
      plt.savefig(options.output + "_efficiency_vs_complexity.png")

  else:
    raise Exception("can not plot complexity scaling")

def plot_efficiency_vs_tolerance (data, options):
  import matplotlib.pyplot as plt

  plt.figure(
      figsize = (
        options.width/100.0,
        options.height/100.0
        ),
      dpi = 100
      )

  tolerance_values = data.get_tolerance()
  if options.thread:
    thread_values = sorted([ int(i) for i in options.thread ])
  else:
    thread_values = data.get_threads()

  if 0 in tolerance_values:
    for t in thread_values:
      walltime = []
      complexity_values = []
      walltime_1 = 0
      for i in range(len(tolerance_values)):
        query = data.get_record(tolerance = tolerance_values[i], threads = t)
        if len(query) == 0:
          raise Exception("can not find result for "
          + "tolerance {:1.3f} and {:d} threads".format(
            tolerance_values[i], t
            )
          )
        walltime.append(query[0]["walltime"])
        complexity_values.append(query[0]["complexity"])
        if tolerance_values[i] == 0:
          walltime_1 = walltime[-1]
      efficiency_values = [ 100*complexity_values[i]*walltime_1/walltime[i] for i in range(len(walltime)) ]
      plt.semilogx(
          tolerance_values,
          efficiency_values,
          linestyle = "-",
          marker = "o",
          label = "{:d} threads".format(t)
          )
      if options.print:
        print("tolerance: ", tolerance_values)
        print("efficiency:", efficiency_values)

    plt.grid(True)
    plt.legend(loc = "lower left")
    plt.xlabel("tolerance")
    plt.ylabel("tolerance effciciency")
    if not options.no_title:
      plt.title("N = {:d}, N_basic = {:d}".format(data.N_chunk, data.N_basic))

    if options.output:
      plt.savefig(options.output + "_efficiency_vs_tolerance.png")

  else:
    raise Exception("can not plot tolerance scaling")

def main ():
  import argparse
  import matplotlib.pyplot as plt

  parser = argparse.ArgumentParser()

  parser.add_argument("FILE",
      help = "The output file of the scaling script",
      nargs = "+")

  parser.add_argument("--no-title",
      help = "Do not add a title to the graphs",
      default = False,
      action = "store_true")

  parser.add_argument("--output",
      help = "Save figures into FILEBASE")

  parser.add_argument("--width",
      help = "The width of the figure in pixels",
      default = 800,
      type = int)

  parser.add_argument("--height",
      help = "The height of the figure in pixels",
      default = 600,
      type = int)

  parser.add_argument("--print",
      help = "Print the data file",
      default = False,
      action = "store_true")

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

  options = parser.parse_args()

  if options.thread:
    options.thread = flatten_list(options.thread)
    print("plotting only threads " + str(options.thread))

  if options.complexity:
    options.complexity = flatten_list(options.complexity)
    print("plotting only complexity " + str(options.complexity))

  for filename in options.FILE:
    data = scaling_data(filename)
    print(data.info())

    if options.print:
      print(str(data))

    plot_walltime_vs_complexity(data, options)
    plot_walltime_vs_tolerance(data, options)
    plot_walltime_vs_threads(data, options)
    plot_walltime(data, options)
    plot_efficiency_vs_threads(data, options)
    plot_efficiency_vs_complexity(data, options)
    plot_efficiency_vs_tolerance(data, options)

  if not options.output:
    plt.show()

if __name__ == "__main__":
  main()
