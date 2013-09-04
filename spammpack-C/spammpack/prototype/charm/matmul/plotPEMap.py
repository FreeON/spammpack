#!/usr/bin/env python

import argparse
import numpy as np
import re
import subprocess
import sys
import tempfile

##############################################

def script (line):
  """ Print a line in the POVRay script. """
  global povray_script
  global python_version
  if python_version.major == 2:
    povray_script.write(line + "\n")
  else:
    povray_script.write(bytes(line + "\n", "UTF-8"))

def openPOVRay ():
  global povray_script
  povray_script = tempfile.NamedTemporaryFile(suffix = ".pov", delete = False)
  print("writing POVRay script into {:s}".format(povray_script.name))

def getColor (value, N):
  return "red 1 green 0 blue 0"

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

def generatePOVRay (iteration, PEMap_A, PEMap_C, PEMap_convolution):
  global povray_script
  openPOVRay()

  ( N, _ ) = PEMap_A.shape

  script("/* Automatically generated... */")
  script("#include \"colors.inc\"")

  # Place the camera.
  script("camera {")
  script("  location  < {:e}, {:e}, {:e} >".format(
    2*N, 2*N, 2*N))
  script("  look_at < 0, 0, 0 >")
  script("}")

  # Plot A map on x-y plane.
  script("box {")
  script("  < 0, 0, -0.2 >, < {:f}, {:f}, -0.2 >".format(
    N, N))
  script("  pigment { color White }")
  script("}")

  for i in range(N):
    for j in range(N):
      script("box {")
      script("  < {:f}, {:f}, 0 >, < {:f}, {:f}, 0 >".format(
        0.1+i, 0.1+j, 0.9+i, 0.9+j))
      color_vector = getColor(PEMap_A[i, j], PEMap_A.size)
      script("  pigment {{ color {:s} }}".format(color_vector))
      script("}")

  # Plot B map on x-z plane.
  script("box {")
  script("  < 0, -0.2, 0 >, < {:f}, -0.2, {:f} >".format(
    N, N))
  script("  pigment { color White }")
  script("}")

  for i in range(N):
    for j in range(N):
      script("box {")
      script("  < {:f}, 0, {:f} >, < {:f}, 0, {:f} >".format(
        0.1+i, 0.1+j, 0.9+i, 0.9+j))
      color_vector = getColor(PEMap_A[i, j], PEMap_A.size)
      script("  pigment {{ color {:s} }}".format(color_vector))
      script("}")

  # Plot C map on y-z plane.
  script("box {")
  script("  < -0.2, 0, 0 >, < -0.2, {:f}, {:f} >".format(
    N, N))
  script("  pigment { color White }")
  script("}")

  for i in range(N):
    for j in range(N):
      script("box {")
      script("  < 0, {:f}, {:f} >, < 0, {:f}, {:f} >".format(
        0.1+i, 0.1+j, 0.9+i, 0.9+j))
      color_vector = getColor(PEMap_C[i, j], PEMap_C.size)
      script("  pigment {{ color {:s} }}".format(color_vector))
      script("}")

  # Plot convolution.
  for i in range(N):
    for j in range(N):
      for k in range(N):
        script("box {")
        script("  < {:f}, {:f}, {:f} >, < {:f}, {:f}, {:f} >".format(
          0.1+i, 0.1+j, 0.1+k, 0.9+i, 0.9+j, 0.9+k))
        color_vector = getColor(
            PEMap_convolution[i, j, k], PEMap_convolution.size)
        script("  pigment {{ color {:s} transmit 0.6 }}".format(color_vector))
        script("  finish { metallic 0.4 }")
        script("  hollow")
        script("}")

  povray_script.close()
  render(iteration, povray_script.name)

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

for line in fd:
  if options.debug:
    print("read: ", line.rstrip())

  result = re.compile("iteration ([0-9]+) on").search(line)
  if result:
    iteration = int(result.group(1))
    print("iteration {:d}".format(iteration))
    continue

  result = re.compile("] PEMap for (.*):").search(line)
  if result:
    mapName = result.group(1)
    if inMap and mapName != currentMap:
      raise(Exception("map {:s} already open for reading".format(mapName)))
    if not inMap:
      N = 0
      elementBuffer = []
    inMap = True
    currentMap = mapName
    if options.debug:
      print("opening map {:s}".format(currentMap))

  result = re.compile("PEMap\(([0-9]+),([0-9]+)\) = ([0-9]+)").search(line)
  if result:
    i = int(result.group(1))
    j = int(result.group(2))
    PE = int(result.group(3))
    if not inMap:
      raise(Exception("no map open for reading"))
    elementBuffer.append( (i, j, PE) )
    if i+1 > N:
      N = i+1
    if j+1 > N:
      N = j+1
    continue

  result = re.compile("PEMap\(([0-9]+),([0-9]+),([0-9]+)\) = ([0-9]+)").search(line)
  if result:
    i = int(result.group(1))
    j = int(result.group(2))
    k = int(result.group(3))
    PE = int(result.group(4))
    if not inMap:
      raise(Exception("no map open for reading"))
    elementBuffer.append( (i, j, k, PE) )
    if i+1 > N:
      N = i+1
    if j+1 > N:
      N = j+1
    if k+1 > N:
      N = k+1
    continue

  result = re.compile("end of PEMap for (.*)").search(line)
  if result:
    if currentMap == "convolution":
      PEMap[currentMap] = np.empty([N, N, N], dtype = np.int16)
      PEMap[currentMap].fill(-1)
      for (i, j, k, PE) in elementBuffer:
        PEMap[currentMap][i,j,k] = PE
      if options.render:
        generatePOVRay(iteration, PEMap["matrix A"], PEMap["matrix C"], PEMap["convolution"])
      if options.printPEMap:
        print(PEMap)
    else:
      PEMap[currentMap] = np.empty([N, N], dtype = np.int16)
      PEMap[currentMap].fill(-1)
      for (i, j, PE) in elementBuffer:
        PEMap[currentMap][i,j] = PE
    inMap = False
    if options.debug:
      print("closing map {:s}".format(currentMap))
    continue

fd.close()
