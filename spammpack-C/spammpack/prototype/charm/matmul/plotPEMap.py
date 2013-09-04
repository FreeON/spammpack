#!/usr/bin/env python

import argparse
import math
import numpy as np
import re
import subprocess
import sys
import tempfile

##############################################

def script (line):
  """ Print a line in the POVRay script. """
  global script_file
  global python_version
  if python_version.major == 2:
    script_file.write(line)
  else:
    script_file.write(bytes(line, "UTF-8"))

def openScriptFile (iteration, suffix):
  global script_file
  script_file = tempfile.NamedTemporaryFile(
      suffix = ".{:d}.{:s}".format(iteration, suffix), delete = False)
  print("writing script into {:s}".format(script_file.name))

def getColor (value, N):
  return [
      math.cos(value/float(N-1)*math.pi/2.),
      math.sin(value/float(N-1)*math.pi/2.),
      0 ]

def render (iteration, filename):
  try:
    cmd = [
        "povray",
        "-d",
        "Verbose=false",
        "+OPEMap_{:d}.png".format(iteration),
        "+H1080",
        "+W1920",
        filename ]
    povray = subprocess.Popen(
        cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    povray.wait()
  except subprocess.CalledProcessError as e:
    print("error spawning povray using: " + e.cmd)

  if povray.returncode != 0:
    print("POVRAY: " + cmd.__str__())
    for line in povray.stdout:
      print("POVRAY: " + line.rstrip().decode())
    for line in povray.stderr:
      print("POVRAY: " + line.rstrip().decode())

def generatePOVRay (iteration, numPEs, PEMap_A, PEMap_C, PEMap_convolution):
  global script_file
  openScriptFile(iteration, "pov")

  ( N, _ ) = PEMap_A.shape

  script("/* Automatically generated... */\n")
  script("#include \"colors.inc\"\n")

  # Place the camera.
  script("camera {\n")
  script("  location  < {:e}, {:e}, {:e} >\n".format(
    2*N, 2*N, 2*N))
  script("  look_at < 0, 0, 0 >\n")
  script("}\n")

  # Plot A map on x-y plane.
  script("box {\n")
  script("  < 0, 0, -0.2 >, < {:f}, {:f}, -0.2 >\n".format(
    N, N))
  script("  pigment { color White }\n")
  script("}\n")

  for i in range(N):
    for j in range(N):
      script("box {\n")
      script("  < {:f}, {:f}, 0 >, < {:f}, {:f}, 0 >\n".format(
        0.1+i, 0.1+j, 0.9+i, 0.9+j))
      color_vector = getColor(PEMap_A[i, j], numPEs)
      script("  pigment {{ color red {:f} green {:f} blue {:f} }}\n".format(
        color_vector[0],
        color_vector[1],
        color_vector[2]))
      script("}\n")

  # Plot B map on x-z plane.
  script("box {\n")
  script("  < 0, -0.2, 0 >, < {:f}, -0.2, {:f} >\n".format(
    N, N))
  script("  pigment { color White }\n")
  script("}\n")

  for i in range(N):
    for j in range(N):
      script("box {\n")
      script("  < {:f}, 0, {:f} >, < {:f}, 0, {:f} >\n".format(
        0.1+i, 0.1+j, 0.9+i, 0.9+j))
      color_vector = getColor(PEMap_A[i, j], numPEs)
      script("  pigment {{ color red {:f} green {:f} blue {:f} }}\n".format(
        color_vector[0],
        color_vector[1],
        color_vector[2]))
      script("}\n")

  # Plot C map on y-z plane.
  script("box {\n")
  script("  < -0.2, 0, 0 >, < -0.2, {:f}, {:f} >\n".format(
    N, N))
  script("  pigment { color White }\n")
  script("}\n")

  for i in range(N):
    for j in range(N):
      script("box {\n")
      script("  < 0, {:f}, {:f} >, < 0, {:f}, {:f} >\n".format(
        0.1+i, 0.1+j, 0.9+i, 0.9+j))
      color_vector = getColor(PEMap_C[i, j], numPEs)
      script("  pigment {{ color red {:f} green {:f} blue {:f} }}\n".format(
        color_vector[0],
        color_vector[1],
        color_vector[2]))
      script("}\n")

  # Plot convolution.
  for i in range(N):
    for j in range(N):
      for k in range(N):
        if PEMap_convolution[i, j, k] >= 0:
          script("box {\n")
          script("  < {:f}, {:f}, {:f} >, < {:f}, {:f}, {:f} >\n".format(
            0.1+i, 0.1+j, 0.1+k, 0.9+i, 0.9+j, 0.9+k))
          color_vector = getColor(
              PEMap_convolution[i, j, k], numPEs)
          script("  pigment {{ color red {:f} green {:f} blue {:f} transmit 0.6 }}\n".format(
            color_vector[0],
            color_vector[1],
            color_vector[2]))
          script("  finish { metallic 0.4 }\n")
          script("  hollow\n")
          script("}\n")

  script_file.close()
  render(iteration, script_file.name)

def generateMathematica (iteration, numPEs, PEMap_A, PEMap_C, PEMap_convolution):
  global script_file
  openScriptFile(iteration, "nb")

  # Plot convolution.
  script("Graphics3D[ {\n");
  for i in range(N):
    for j in range(N):
      for k in range(N):
        if PEMap_convolution[i, j, k] >= 0:
          color_vector = getColor(
              PEMap_convolution[i, j, k], numPEs)
          script("  RGBColor[ {:f}, {:f}, {:f} ], ".format(
                color_vector[0],
                color_vector[1],
                color_vector[2]))
          script("Opacity[ 0.1 ], Cuboid[")
          script("{{ {:f}, {:f}, {:f} }}, ".format(0.1+i, 0.1+j, 0.1+k))
          script("{{ {:f}, {:f}, {:f} }} ],\n".format(0.9+i, 0.9+j, 0.9+k))
  script("} ]\n")

  script_file.close()

##############################################

global python_version
python_version = sys.version_info

if python_version.major == 2 and python_version.minor < 7:
  print("need at least python 2.7 (running {:d}.{:d})".format(
    python_version.major, python_version.minor))
  sys.exit(1)

parser = argparse.ArgumentParser()

parser.add_argument("FILE",
    help = "The output file to plot. A value of '-' means standard input.")

parser.add_argument("--render",
    help = "Render the PEMaps",
    action = "store_true",
    default = False)

parser.add_argument("--print",
    help = "Print the PEMaps to stdout",
    dest = "printPEMap",
    action = "store_true",
    default = False)

parser.add_argument("--mathematica",
    help = "Generate Mathematic statements",
    action = "store_true",
    default = False)

parser.add_argument("--tolerance",
    help = "When printing the convolution, filter with TOLERANCE",
    type = float,
    default = 0)

parser.add_argument("--debug",
    help = "Print debugging stuff",
    action = "store_true",
    default = False)

options = parser.parse_args()

if options.FILE == "-":
  fd = sys.stdin
else:
  fd = open(options.FILE)

iteration = -1

inMap = False
currentMap = ""

PEMap = {}
numPEs = -1
line_number = 0

for line in fd:
  line_number += 1
  if options.debug:
    print("read ({:d})".format(line_number), line.rstrip())

  result = re.compile("iteration ([0-9]+) on").search(line)
  if result:
    iteration = int(result.group(1))
    print("iteration {:d}".format(iteration))
    continue

  result = re.compile("] PEMap for (.*):").search(line)
  if result:
    mapName = result.group(1)
    if inMap and mapName != currentMap:
      raise(Exception("line {:d}: map {:s} already open for reading".format(
        line_number, mapName)))
    if not inMap:
      N = 0
      elementBuffer = []
    inMap = True
    currentMap = mapName
    if options.debug:
      print("opening map {:s}".format(currentMap))

  result = re.compile("PEMap\(([0-9]+),([0-9]+)\) = ([0-9]+) \(norm = ([0-9.e+-]+)\)").search(line)
  if result:
    i = int(result.group(1))
    j = int(result.group(2))
    PE = int(result.group(3))
    norm = float(result.group(4))
    if not inMap:
      raise(Exception("line {:d}: no map open for reading".format(line_number)))
    elementBuffer.append( (i, j, PE, norm) )
    if i+1 > N:
      N = i+1
    if j+1 > N:
      N = j+1
    if PE+1 > numPEs:
      numPEs = PE+1
    continue

  result = re.compile("PEMap\(([0-9]+),([0-9]+),([0-9]+)\) = ([0-9]+) \(norm = ([0-9.e+-]+)\)").search(line)
  if result:
    i = int(result.group(1))
    j = int(result.group(2))
    k = int(result.group(3))
    PE = int(result.group(4))
    norm = float(result.group(5))
    if not inMap:
      raise(Exception("line {:d}: no map open for reading".format(line_number)))
    elementBuffer.append( (i, j, k, PE, norm) )
    if i+1 > N:
      N = i+1
    if j+1 > N:
      N = j+1
    if k+1 > N:
      N = k+1
    if PE+1 > numPEs:
      numPEs = PE+1
    continue

  result = re.compile("end of PEMap for (.*)").search(line)
  if result:
    if currentMap == "convolution":
      PEMap[currentMap] = np.empty([N, N, N], dtype = np.int16)
      PEMap[currentMap].fill(-1)
      for (i, j, k, PE, norm) in elementBuffer:
        if norm > options.tolerance:
          PEMap[currentMap][i,j,k] = PE
      if options.printPEMap:
        print(PEMap)
      if options.render:
        generatePOVRay(
            iteration, numPEs, PEMap["matrix A"], PEMap["matrix C"],
            PEMap["convolution"])
      if options.mathematica:
        generateMathematica(
            iteration, numPEs, PEMap["matrix A"], PEMap["matrix C"],
            PEMap["convolution"])
    else:
      PEMap[currentMap] = np.empty([N, N], dtype = np.int16)
      PEMap[currentMap].fill(-1)
      for (i, j, PE, norm) in elementBuffer:
        PEMap[currentMap][i,j] = PE
    inMap = False
    if options.debug:
      print("closing map {:s}".format(currentMap))
    continue

fd.close()
