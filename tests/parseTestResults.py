#!/usr/bin/python
#
# This code is part of the MondoSCF suite of programs for linear scaling
# electronic structure theory and ab initio molecular dynamics.
#
# Copyright (2004). The Regents of the University of California. This
# material was produced under U.S. Government contract W-7405-ENG-36
# for Los Alamos National Laboratory, which is operated by the University
# of California for the U.S. Department of Energy. The U.S. Government has
# rights to use, reproduce, and distribute this software.  NEITHER THE
# GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
# OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version. Accordingly, this program is distributed in
# the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License at www.gnu.org for details.
#
# While you may do as you like with this software, the GNU license requires
# that you clearly mark derivative software.  In addition, you are encouraged
# to return derivative works to the MondoSCF group for review, and possible
# disemination in future releases.
#
# Regression tests for FreeON.
#
# Nicolas Bock <nbock@lanl.gov>

import argparse
import logging
import math
import os
import re
import sys

########################################
#
# Class definitions.
#
########################################

class referenceTag:

  KnownTagTypes = [ "equal", "less" ]

  def __init__ (self, tagType):
    if tagType in referenceTag.KnownTagTypes:
      self.tagType = tagType
    else:
      log.error("unknown tag type (%s)" % (tagType))
      sys.exit(1)
    self.tag = None
    self.values = []
    self.valueFromOutput = []
    self.foundMatch = []
    self.result = []
    self.unreferencedMatches = []
    self.absoluteTolerance = 0.0
    self.relativeTolerance = 0.0
    self.numberFailed = 0

  def setTag (self, tagValue):
    self.tag = tagValue

  def setAbsoluteTolerance (self, absoluteTolerance):
    self.absoluteTolerance = math.fabs(absoluteTolerance)
    self.relativeTolerance = 0.0

  def setRelativeTolerance (self, relativeTolerance):
    self.absoluteTolerance = 0.0
    self.relativeTolerance = math.fabs(relativeTolerance)

  def addValue (self, value):
    self.values.append(value)
    self.valueFromOutput.append(False)
    self.foundMatch.append(False)
    self.result.append(False)

  def check (self, line):
    result = re.compile(self.tag).search(line)
    if result:
      convertedValue = re.sub("[dD]", "e", result.group(1))
      try:
        value = float(convertedValue)
      except ValueError:
        log.error("illegal value, regular expression \"%s\", line = \"%s\", value = \"%s\""
            % (self.tag, line, convertedValue))
        sys.exit(1)

      if self.tagType == "equal":
        haveReference = False
        for i in range(len(self.foundMatch)):
          if not self.foundMatch[i]:
            haveReference = True
            self.foundMatch[i] = True
            self.valueFromOutput[i] = value
            if self.absoluteTolerance > 0.0:
              if math.fabs(value-self.values[i]) <= self.absoluteTolerance:
                self.result[i] = True
              else:
                self.numberFailed += 1
            elif self.relativeTolerance > 0.0:
              if math.fabs((value-self.values[i])/self.values[i]) <= self.relativeTolerance:
                self.result[i] = True
              else:
                self.numberFailed += 1
            else:
              if value == self.values[i]:
                self.result[i] = True
              else:
                self.numberFailed += 1
            break
        if not haveReference:
          self.unreferencedMatches.append(value)

      elif self.tagType == "less":
        haveReference = False
        for i in range(len(self.foundMatch)):
          if not self.foundMatch[i]:
            haveReference = True
            self.foundMatch[i] = True
            self.valueFromOutput[i] = value
            if self.absoluteTolerance > 0.0:
              if value <= self.values[i]+self.absoluteTolerance:
                self.result[i] = True
              else:
                self.numberFailed += 1
            elif self.relativeTolerance > 0.0:
              if value <= (1+self.relativeTolerance)*self.values[i]:
                self.result[i] = True
              else:
                self.numberFailed += 1
            else:
              if value <= self.values[i]:
                self.result[i] = True
              else:
                self.numberFailed += 1
            break
        if not haveReference:
          self.unreferencedMatches.append(value)

  def printResult (self):
    resultString = "type: %s\n" % (self.tagType)
    resultString += "regular expression: \"%s\"\n" % (self.tag)
    if self.relativeTolerance > 0.0:
      resultString += "relative tolerance: %1.2e\n" % (self.relativeTolerance)
    else:
      resultString += "absolute tolerance: %e\n" % (self.absoluteTolerance)
    for i in range(len(self.values)):
      if self.foundMatch[i]:
        resultString += "(value found: %e)" % (self.valueFromOutput[i])
        if self.tagType == "equal":
          if self.result[i]:
            resultString += " == "
          else:
            resultString += " != "
        elif self.tagType == "less":
          if self.result[i]:
            resultString += " <= "
          else:
            resultString += " >  "
        resultString += "(reference: %e)" % (self.values[i])
        if self.result[i]:
          resultString += "\n"
        else:
          resultString += ", outside of tolerance, "
          if self.relativeTolerance > 0.0:
            resultString += "|diff/ref| = %1.4e\n" % (math.fabs((self.values[i]
              -self.valueFromOutput[i])/self.values[i]))
          else:
            resultString += "|diff| = %e\n" % (math.fabs(self.values[i]-self.valueFromOutput[i]))
      else:
        resultString += "did not find a matching tag in output for reference %e\n" % (self.values[i])
    resultString += "unreferenced matches: %d\n" % (len(self.unreferencedMatches))
    for match in self.unreferencedMatches:
      resultString += "unreferenced match value: %e\n" % (match)
    return resultString

  def __repr__ (self):
    return "[%s, %d values, abs tol %1.2e, rel tol %1.2e]" % (self.tagType,
        len(self.values), self.absoluteTolerance, self.relativeTolerance)

  def __str__ (self):
    stringResult = ""

    stringResult += self.tagType + "\n{\n"
    stringResult += "  " + "tag = \"%s\";\n" % (self.tag)
    if self.absoluteTolerance > 0.0:
      stringResult += "  " + "absoluteTolerance = %1.12e;\n" % (self.absoluteTolerance)
    if self.relativeTolerance > 0.0:
      stringResult += "  " + "relativeTolerance = %1.12e;\n" % (self.relativeTolerance)
    for i in self.valueFromOutput:
      stringResult += "  " + "value = %1.12e;\n" % (i)
    stringResult += "}"

    return stringResult

########################################
#
# Main program.
#
########################################

parser = argparse.ArgumentParser()

parser.add_argument("--reference",
    metavar = "FILE",
    help = "load reference values from FILE",
    dest = "reference",
    required = True)

parser.add_argument("--output",
    metavar = "FILE",
    help = "load FreeON output from FILE instead of from standard input",
    dest = "output")

parser.add_argument("--verbose", "-v",
    help = "print lots of stuff.",
    dest = "verbose",
    action = "store_true",
    default = False)

parser.add_argument("--log-output",
    help = "log the output to the log file",
    action = "store_true",
    default = False)

parser.add_argument("--log-filename",
    metavar = "FILE",
    help = "log to FILE (default = %(default)s)",
    default = "regressionTest.log")

options = parser.parse_args()

# Set up logging.
log = logging.getLogger("regressionTest")
log.setLevel(logging.DEBUG)

# Set console logger.
logHandler = logging.StreamHandler()
logFormatter = logging.Formatter("[%(name)s] %(message)s")
logHandler.setFormatter(logFormatter)
logHandler.setLevel(logging.INFO)

if options.verbose:
  logHandler.setLevel(logging.DEBUG)

log.addHandler(logHandler)

# Set file logger.
logHandler = logging.FileHandler(options.log_filename)
logFormatter = logging.Formatter("%(asctime)s [%(name)s %(levelname)s] %(message)s", "%y-%m-%d %H:%M:%S")
logHandler.setFormatter(logFormatter)
logHandler.setLevel(logging.DEBUG)

log.addHandler(logHandler)

log.debug("starting new test...")

if not options.output:
  # Load output from standard input.
  log.debug("loading output from standard input")
  output = sys.stdin.readlines()

else:
  # Load output from file.
  log.debug("loading output from file: " + options.output)
  fd = open(options.output)
  output = fd.readlines()
  fd.close()

# Parse output against reference.
log.debug("loading reference tags")

fd = open(options.reference)
lines = fd.readlines()
fd.close()

referenceString = ""
for i in range(len(lines)):
  line = lines[i]

  # Remove comments.
  tokens = line.split('#')
  line = tokens[0].strip()
  if len(line) == 0:
    continue

  referenceString += line

# Interpret the reference.
reference = []
i = 0
while i < len(referenceString):
  log.debug("i = %d" % (i))

  substring = referenceString[i:len(referenceString)]
  log.debug("substring = %s" % (substring))

  blockStart = substring.find("{")
  if blockStart < 0:
    log.error("missing opening block")
    sys.exit(1)

  blockEnd = substring.find("}")
  if blockEnd < 0:
    log.error("missing closing block")
    sys.exit(1)

  tagName = substring[0:blockStart].strip()
  log.debug("found tag %s" % (tagName))

  blockString = substring[blockStart+1:blockEnd]
  log.debug("block %s" % (blockString))

  newTagObject = referenceTag(tagName)

  tokens = blockString.split(";")
  if len(tokens[len(tokens)-1]) != 0:
    log.error("syntax error: missing \";\" at the end of block")
    sys.exit(1)

  for token in tokens:
    if len(token) == 0:
      continue

    if token.find("=") < 0:
      log.error("syntax error (missing =): %s" % (token))
      sys.exit(1)

    key = token[0:token.find("=")-1].strip()
    value = token[token.find("=")+1:len(token)].strip()
    log.debug("key = %s" % (key))
    log.debug("value = %s" % (value))

    if key == "tag":
      if value[0] != "\"" or value[len(value)-1] != "\"":
        log.error("tag regular expressions need to be enclosed in \"")
        sys.exit(1)
      value = value[1:len(value)-1]
      newTagObject.setTag(value)

    elif key == "value":
      value = re.sub("[dD]", "e", value)
      try:
        newTagObject.addValue(float(value))
      except ValueError:
        log.error("illegal number: \"%s\"" % (value))
        sys.exit(1)

    elif key == "absoluteTolerance":
      value = re.sub("[dD]", "e", value)
      try:
        newTagObject.setAbsoluteTolerance(float(value))
      except ValueError:
        log.error("illegal number: \"%s\"" % (value))
        sys.exit(1)

    elif key == "relativeTolerance":
      value = re.sub("[dD]", "e", value)
      try:
        newTagObject.setRelativeTolerance(float(value))
      except ValueError:
        log.error("illegal number: \"%s\"" % (value))
        sys.exit(1)

    else:
      log.error("syntax error (unknown key): %s" % (token))
      sys.exit(1)

  reference.append(newTagObject)

  i += blockEnd+1

log.debug("checking tags: " + str(reference))

successfullyTerminated = False
scratchDirectory = "unknown"

linenumber = 0
inputFile = None
for line in output:
  linenumber += 1
  line = line.rstrip()

  if options.log_output:
    log.debug(line.rstrip())

  check = re.compile("Successful FreeON run").search(line)
  if check:
    successfullyTerminated = True

  check = re.compile("scratch directory at (.*)").search(line)
  if check:
    scratchDirectory = check.group(1)

  check = re.compile("input file = (.*)").search(line)
  if check:
    inputFile = check.group(1)

  for tag in reference:
    tag.check(line)

exitStatus = 0

if inputFile:
  log.debug("input file: " + inputFile)
if options.output:
  log.debug("output file: " + options.output)
else:
  log.debug("reading from standard input")
log.debug("reference file: " + options.reference)
log.debug("scratch directory: " + scratchDirectory)
log.debug("log file: " + os.path.join(os.getcwd(), "regressionTest.log"))
if successfullyTerminated:
  log.debug("FreeON successfully terminated")
else:
  exitStatus = 1
  log.debug("FreeON did not successfully terminate")

log.debug("result:")
numberFailed = 0
numberUnreferencedMatches = 0
for tag in reference:
  numberFailed += tag.numberFailed
  numberUnreferencedMatches += len(tag.unreferencedMatches)
  resultString = tag.printResult()
  tokens = resultString.split("\n")
  for token in tokens:
    log.debug("  " + token)
  log.debug("")

log.debug("%d failed match(es)" % (numberFailed))
log.debug("%d unreferenced match(es)" % (numberUnreferencedMatches))

log.debug("")

if numberFailed > 0:
  log.debug("if your reference file had looked like this, I would")
  log.debug("have been happier...")
  # Print reference file with proper values.
  for i in reference:
    log.debug("\n" + i.__str__())
  log.debug("")

# Exit with a proper exit code.
if numberFailed == 0:
  log.debug("successful test")
else:
  log.debug("test failed")
  exitStatus = 1

sys.exit(exitStatus)
