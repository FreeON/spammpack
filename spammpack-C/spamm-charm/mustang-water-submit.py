#!/usr/bin/env python

def string_compat (line):
  import sys
  if sys.version_info.major > 2:
    return bytes(line, encoding = "UTF-8")
  else:
    return line

def main ():
  import argparse
  import logging
  import re
  import subprocess
  import tempfile

  parser = argparse.ArgumentParser()

  parser.add_argument(
      "--spamm-version",
      help = "The SpAMM version (default is %(default)s)",
      choices = [ "serial", "openmp" ],
      default = "openmp"
      )

  parser.add_argument(
      "--name",
      help = "The name of the water input files",
      required = True
      )

  parser.add_argument(
      "--fockians",
      help = "The path to the Fockians (default is %(default)s)",
      default = "/scratch/nbock"
      )

  parser.add_argument(
      "--template",
      help = "The location of the job template (default is %(default)s)",
      default = "mustang.job.water.template"
      )

  parser.add_argument(
      "--Ne",
      help = "The number of electrons for SP2",
      type = int,
      required = True
      )

  parser.add_argument(
      "--tolerance",
      help = "The SpAMM tolerance (default is %(default)s)",
      type = float,
      default = 0
      )

  parser.add_argument(
      "--blocksize",
      help = "The blocksize of the chunk (default is %(default)s)",
      type = int,
      default = 128
      )

  parser.add_argument(
      "--basic",
      help = "The size of the basic submatrices (default is %(default)s)",
      type = int,
      default = 4
      )

  parser.add_argument(
      "--nodes",
      help = "The number of nodes to allocate (default is %(default)s)",
      type = int,
      default = 1
      )

  parser.add_argument(
      "--threads",
      help = "The number of threads (default is not to restrict)",
      type = int,
      default = -1
      )

  parser.add_argument(
      "--no-submit",
      help = "do not submit job",
      default = False,
      action = "store_true"
      )

  options = parser.parse_args()

  logging.basicConfig(level = logging.INFO)

  template = open(options.template)
  jobscript = tempfile.NamedTemporaryFile(
      prefix = "mustang-water.",
      suffix = ".job",
      dir = ".",
      delete = False)

  logging.info("writing into {:s}".format(jobscript.name))

  for line in template:
    line = line.rstrip()

    # Parse and replace.
    line = re.sub("nodes=N", "nodes=%d" % (options.nodes), line)
    if options.threads > 0:
      line = re.sub("OMP_NUM_THREADS=1", "OMP_NUM_THREADS=%d" % (options.threads), line)
    else:
      line = re.sub("export OMP", "#export OMP", line)
    line = re.sub("JOBNAME", "%s-%s" % (options.name, options.spamm_version), line)
    line = re.sub("SPAMM_VERSION", options.spamm_version, line)
    line = re.sub("NE", "%d" % (options.Ne), line)
    line = re.sub("FOCKIAN", options.fockians + "/" + options.name + ".F_DIIS", line)
    line = re.sub("TOLERANCE", "%e" % (options.tolerance), line)
    line = re.sub("BLOCK", "%d" % (options.blocksize), line)
    line = re.sub("BASIC", "%d" % (options.basic), line)

    # Write new line.
    jobscript.write(string_compat(line + "\n"))

  template.close()
  jobscript.close()

  if options.no_submit:
    logging.info("not submitting job")
  else:
    logging.info("submitting job")
    msub = subprocess.Popen([ "msub", jobscript.name ])
    msub.wait()

if __name__ == "__main__":
  main()
