#!/usr/bin/python
#
# vim: tw=0

"""Generate SSE assembly code for a kernel operating on a 4x4 blocks.

The script generates a kernel written in assembly using SSE instructions that
operates on basic 4x4 matrix blocks. The kernel tier can be specified in the
function call.
"""

import logging
import math
import optparse
import os.path
import sys

class SSERegister:
  """A class that takes care of returning an available SSE register to the
  caller. This is how we rotate through the 16 available registers on an x86_64
  CPU."""

  variables = {}
  registerPool = [ "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm4", "%xmm5",
      "%xmm6", "%xmm7", "%xmm8", "%xmm9", "%xmm10", "%xmm11", "%xmm12",
      "%xmm13", "%xmm14", "%xmm15" ]

  def __init__ (self, name, registerName = None):
    """Get a free register for the named variable. The register needs to be
    released again to become free again."""

    if name in SSERegister.variables:
      log.error("The variable %s is already assigned a register" % (name))
      sys.exit(1)

    if registerName and not registerName in SSERegister.registerPool:
      log.error("The requested register %s is not available anymore" % (registerName))
      sys.exit(1)

    if registerName:
      self.name = name
      self.register = registerName
      SSERegister.registerPool.remove(registerName)
      SSERegister.variables[name] = registerName
    else:
      if len(SSERegister.registerPool) > 0:
        self.name = name
        self.register = SSERegister.registerPool.pop(0)
        SSERegister.variables[name] = self.register
      else:
        log.error("no registers left")
        sys.exit(1)

    log.debug("assigned %s --> %s" % (self.register, self.name))
    log.debug("registerPool: %s" % (SSERegister.registerPool))
    log.debug("variables: %s" % (SSERegister.variables))

  def release (self):
    """Release a register back to the register pool."""

    if not self.name in SSERegister.variables:
      log.error("The variable %s is not assigned a register" % (self.name))
      sys.exit(1)

    SSERegister.registerPool.append(self.register)
    del SSERegister.variables[self.name]

    log.debug("released register %s assigned to variable %s" % (self.register, self.name))
    log.debug("registerPool: %s" % (SSERegister.registerPool))
    log.debug("variables: %s" % (SSERegister.variables))

  def getName (self):
    """Return the name of the variable."""
    return self.name

  def getRegister (self):
    """Return the register assigned to this variable."""
    return self.register

  def __str__ (self):
    """Use this variable in the code."""

    if not self.name in SSERegister.variables:
      log.error("The variable %s is not assigned a register" % (self.name))
      sys.exit(1)
    return self.register

class counter:
  """A counter object."""

  def __init__ (self):
    """Initialize the counter to zero."""
    self.counter = 0

  def __init__ (self, initial_value):
    """Initialize the counter and reset it to an initial value."""
    self.counter = initial_value

  def increment (self):
    """Increment the counter by one."""
    self.counter += 1

  def get (self):
    """Get the current value of the counter."""
    return self.counter

def Z_curve_index (i, j):
  """Return the linear index of a Z-curve ordered matrix."""

  result = 0
  index = 0
  while i >= 2**index or j >= 2**index:
    if j & 2**index != 0:
      result += 2**(2*index)

    if i & 2**index != 0:
      result += 2**(2*index+1)

    index += 1

  return result

def row_major_index (i, j, N):
  """Return the index within a matrix block."""

  # Row-major storage.
  return i*N+j

def block_product (i, k, j, clearC = True, writeC = True):
  # Fix up indices.
  i -= 1
  j -= 1
  k -= 1

  # Reference the correct C registers.
  global C1
  global C2
  global C3
  global C4

  if clearC:
    C1 = SSERegister("C1")
    C2 = SSERegister("C2")
    C3 = SSERegister("C3")
    C4 = SSERegister("C4")

    print("")
    print("  # Reset C(%d,%d) matrix block accumulators." % (i+1, j+1))
    print("  xorps %s, %s" % (C1, C1))
    print("  xorps %s, %s" % (C2, C2))
    print("  xorps %s, %s" % (C3, C3))
    print("  xorps %s, %s" % (C4, C4))

  if options.generate_checks:
    print("")
    print("  .balign 16")
    print("jump_%d:" % (block_counter.get()))
    block_counter.increment()

    norm = SSERegister("norm")

    print("  # Check norm of product ||A(%d,%d)||*||B(%d,%d)||." % (i+1, k+1, k+1, j+1))
    print("  movss 0x%x(A), %s" % ((i*options.N+k)*4+offset_norm, norm))
    print("  mulss 0x%x(B), %s" % ((k*options.N+j)*4+offset_norm, norm))

    # When comparing with the Intel Software Developer's Manual, keep in
    # mind that Intel uses Intel syntax and this code is writting using
    # At&T syntax, which means that the order of operand 1 and 2 are the
    # opposite.
    print("  comiss %s, %s" % (tolerance, norm))
    print("  jbe jump_%d" % (block_counter.get()))

    norm.release()

  if options.SSE == 1:
    B1 = SSERegister("B1")
    B2 = SSERegister("B2")
    B3 = SSERegister("B3")
    B4 = SSERegister("B4")

    print("")
    print("  # Calculate C(%d,%d) += A(%d,%d)*B(%d,%d)." % (i+1, j+1, i+1, k+1, k+1, j+1))
    print("  movaps 0x%x(B), %s" % (row_major_index(0, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense, B1))
    print("  movaps 0x%x(B), %s" % (row_major_index(1, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense, B2))
    print("  movaps 0x%x(B), %s" % (row_major_index(2, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense, B3))
    print("  movaps 0x%x(B), %s" % (row_major_index(3, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense, B4))

    A11 = SSERegister("A11")
    A12 = SSERegister("A12")
    A13 = SSERegister("A13")

    print("  movaps 0x%x(A), %s" % (row_major_index(0, 0, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A11))
    print("  movaps 0x%x(A), %s" % (row_major_index(0, 1, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A12))
    print("  movaps 0x%x(A), %s" % (row_major_index(0, 2, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A13))
    print("  mulps %s, %s" % (B1, A11))
    print("  mulps %s, %s" % (B2, A12))
    print("  addps %s, %s" % (A11, C1))

    A11.release()
    A14 = SSERegister("A14")

    print("  movaps 0x%x(A), %s" % (row_major_index(0, 3, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A14))
    print("  mulps %s, %s" % (B3, A13))
    print("  addps %s, %s" % (A12, C1))

    A12.release()
    A21 = SSERegister("A21")

    print("  movaps 0x%x(A), %s" % (row_major_index(1, 0, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A21))
    print("  mulps %s, %s" % (B4, A14))
    print("  addps %s, %s" % (A13, C1))

    A13.release()
    A22 = SSERegister("A22")

    print("  movaps 0x%x(A), %s" % (row_major_index(1, 1, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A22))
    print("  mulps %s, %s" % (B1, A21))
    print("  addps %s, %s" % (A14, C1))

    A14.release()
    A23 = SSERegister("A23")

    print("  movaps 0x%x(A), %s" % (row_major_index(1, 2, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A23))
    print("  mulps %s, %s" % (B2, A22))
    print("  addps %s, %s" % (A21, C2))

    A21.release()
    A24 = SSERegister("A24")

    print("  movaps 0x%x(A), %s" % (row_major_index(1, 3, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A24))
    print("  mulps %s, %s" % (B3, A23))
    print("  addps %s, %s" % (A22, C2))

    A22.release()
    A31 = SSERegister("A31")

    print("  movaps 0x%x(A), %s" % (row_major_index(2, 0, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A31))
    print("  mulps %s, %s" % (B4, A24))
    print("  addps %s, %s" % (A23, C2))

    A23.release()
    A32 = SSERegister("A32")

    print("  movaps 0x%x(A), %s" % (row_major_index(2, 1, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A32))
    print("  mulps %s, %s" % (B1, A31))
    print("  addps %s, %s" % (A24, C2))

    A24.release()
    A33 = SSERegister("A33")

    print("  movaps 0x%x(A), %s" % (row_major_index(2, 2, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A33))
    print("  mulps %s, %s" % (B2, A32))
    print("  addps %s, %s" % (A31, C3))

    A31.release()
    A34 = SSERegister("A34")

    print("  movaps 0x%x(A), %s" % (row_major_index(2, 3, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A34))
    print("  mulps %s, %s" % (B3, A33))
    print("  addps %s, %s" % (A32, C3))

    A32.release()
    A41 = SSERegister("A41")

    print("  movaps 0x%x(A), %s" % (row_major_index(3, 0, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A41))
    print("  mulps %s, %s" % (B4, A34))
    print("  addps %s, %s" % (A33, C3))

    A33.release()
    A42 = SSERegister("A42")

    print("  movaps 0x%x(A), %s" % (row_major_index(3, 1, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A42))
    print("  mulps %s, %s" % (B1, A41))
    print("  addps %s, %s" % (A34, C3))

    A34.release()
    A43 = SSERegister("A43")

    print("  movaps 0x%x(A), %s" % (row_major_index(3, 2, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A43))
    print("  mulps %s, %s" % (B2, A42))
    print("  addps %s, %s" % (A41, C4))

    A41.release()
    A44 = SSERegister("A44")

    print("  movaps 0x%x(A), %s" % (row_major_index(3, 3, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A44))
    print("  mulps %s, %s" % (B3, A43))
    print("  addps %s, %s" % (A42, C4))
    print("  mulps %s, %s" % (B4, A44))
    print("  addps %s, %s" % (A43, C4))
    print("  addps %s, %s" % (A44, C4))

    A42.release()
    A43.release()
    A44.release()

    B1.release()
    B2.release()
    B3.release()
    B4.release()

  elif options.SSE == 3:
    print("")

  elif options.SSE == 4.1:
    print("")
    print("  # Calculate C(%d,%d) += A(%d,%d)*B(%d,%d)." % (i+1, j+1, i+1, k+1, k+1, j+1))

    A1 = SSERegister("A1")
    print("  movaps 0x%x(A), %s" % (row_major_index(0, 0, 4)*4+row_major_index(i, k, options.N)*16*4+offset_block_dense, A1))

    B1 = SSERegister("B1")
    print("  movaps 0x%x(B), %s" % (row_major_index(0, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense_transpose, B1))
    B2 = SSERegister("B2")
    print("  movaps 0x%x(B), %s" % (row_major_index(1, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense_transpose, B2))
    B3 = SSERegister("B3")
    print("  movaps 0x%x(B), %s" % (row_major_index(2, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense_transpose, B3))
    B4 = SSERegister("B4")
    print("  movaps 0x%x(B), %s" % (row_major_index(3, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense_transpose, B4))

    print("")
    print("  # Calculate C(1,:).")

    C11 = SSERegister("C11")
    print("  movaps %s, %s" % (B1, C11))
    print("  dpps $0xf1, %s, %s" % (A1, C11))

    C12 = SSERegister("C12")
    print("  movaps %s, %s" % (B2, C12))
    print("  dpps $0xf2, %s, %s" % (A1, C12))

    C13 = SSERegister("C13")
    print("  movaps %s, %s" % (B3, C13))
    print("  dpps $0xf4, %s, %s" % (A1, C13))

    C14 = SSERegister("C14")
    print("  movaps %s, %s" % (B4, C14))
    print("  dpps $0xf8, %s, %s" % (A1, C14))

    print("  blendps $0x01, %s, %s" % (C11, C12))
    C11.release()
    print("  blendps $0x03, %s, %s" % (C12, C13))
    C12.release()
    print("  blendps $0x07, %s, %s" % (C13, C14))
    C13.release()
    print("  addps %s, %s" % (C14, C1))
    C14.release()
    A1.release()

    A2 = SSERegister("A2")
    print("  movaps 0x%x(A), %s" % (row_major_index(1, 0, 4)*4+row_major_index(i, k, options.N)*16*4+offset_block_dense, A2))

    print("")
    print("  # Calculate C(2,:).")

    C21 = SSERegister("C21")
    print("  movaps %s, %s" % (B1, C21))
    print("  dpps $0xf1, %s, %s" % (A2, C21))

    C22 = SSERegister("C22")
    print("  movaps %s, %s" % (B2, C22))
    print("  dpps $0xf2, %s, %s" % (A2, C22))

    C23 = SSERegister("C23")
    print("  movaps %s, %s" % (B3, C23))
    print("  dpps $0xf4, %s, %s" % (A2, C23))

    C24 = SSERegister("C24")
    print("  movaps %s, %s" % (B4, C24))
    print("  dpps $0xf8, %s, %s" % (A2, C24))

    print("  blendps $0x01, %s, %s" % (C21, C22))
    C21.release()
    print("  blendps $0x03, %s, %s" % (C22, C23))
    C22.release()
    print("  blendps $0x07, %s, %s" % (C23, C24))
    C23.release()
    print("  addps %s, %s" % (C24, C2))
    C24.release()
    A2.release()

    A3 = SSERegister("A3")
    print("  movaps 0x%x(A), %s" % (row_major_index(2, 0, 4)*4+row_major_index(i, k, options.N)*16*4+offset_block_dense, A3))

    print("")
    print("  # Calculate C(3,:).")

    C31 = SSERegister("C31")
    print("  movaps %s, %s" % (B1, C31))
    print("  dpps $0xf1, %s, %s" % (A3, C31))

    C32 = SSERegister("C32")
    print("  movaps %s, %s" % (B2, C32))
    print("  dpps $0xf2, %s, %s" % (A3, C32))

    C33 = SSERegister("C33")
    print("  movaps %s, %s" % (B3, C33))
    print("  dpps $0xf4, %s, %s" % (A3, C33))

    C34 = SSERegister("C34")
    print("  movaps %s, %s" % (B4, C34))
    print("  dpps $0xf8, %s, %s" % (A3, C34))

    print("  blendps $0x01, %s, %s" % (C31, C32))
    C31.release()
    print("  blendps $0x03, %s, %s" % (C32, C33))
    C32.release()
    print("  blendps $0x07, %s, %s" % (C33, C34))
    C33.release()
    print("  addps %s, %s" % (C34, C3))
    C34.release()
    A3.release()

    A4 = SSERegister("A4")
    print("  movaps 0x%x(A), %s" % (row_major_index(3, 0, 4)*4+row_major_index(i, k, options.N)*16*4+offset_block_dense, A4))

    print("")
    print("  # Calculate C(4,:).")

    C41 = SSERegister("C41")
    print("  movaps %s, %s" % (B1, C41))
    print("  dpps $0xf1, %s, %s" % (A4, C41))

    C42 = SSERegister("C42")
    print("  movaps %s, %s" % (B2, C42))
    print("  dpps $0xf2, %s, %s" % (A4, C42))

    C43 = SSERegister("C43")
    print("  movaps %s, %s" % (B3, C43))
    print("  dpps $0xf4, %s, %s" % (A4, C43))

    C44 = SSERegister("C44")
    print("  movaps %s, %s" % (B4, C44))
    print("  dpps $0xf8, %s, %s" % (A4, C44))

    print("  blendps $0x01, %s, %s" % (C41, C42))
    C41.release()
    print("  blendps $0x03, %s, %s" % (C42, C43))
    C42.release()
    print("  blendps $0x07, %s, %s" % (C43, C44))
    C43.release()
    print("  addps %s, %s" % (C44, C4))
    C44.release()
    A4.release()

    # Done.
    B1.release()
    B2.release()
    B3.release()
    B4.release()

  if writeC:
    if options.generate_checks:
      print("")
      print("  .balign 16")
      print("jump_%d:" % (block_counter.get()))
      block_counter.increment()

    if not options.alphaOne:
      print("")
      print("  # Multiply C(%d,%d) by alpha." % (i+1, j+1))
      print("  mulps alpha, %s" % (C1))
      print("  mulps alpha, %s" % (C2))
      print("  mulps alpha, %s" % (C3))
      print("  mulps alpha, %s" % (C4))

    print("")
    print("  # Add accumulated C(%d,%d) to already existing." % (i+1, j+1))
    print("  addps 0x%x(C), %s" % (row_major_index(0, 0, 4)*4+row_major_index(i, j, options.N)*16*4+offset_block_dense, C1))
    print("  addps 0x%x(C), %s" % (row_major_index(1, 0, 4)*4+row_major_index(i, j, options.N)*16*4+offset_block_dense, C2))
    print("  addps 0x%x(C), %s" % (row_major_index(2, 0, 4)*4+row_major_index(i, j, options.N)*16*4+offset_block_dense, C3))
    print("  addps 0x%x(C), %s" % (row_major_index(3, 0, 4)*4+row_major_index(i, j, options.N)*16*4+offset_block_dense, C4))

    print("")
    print("  # Write out C(%d,%d) submatrix block." % (i+1, j+1))
    print("  movaps %s, 0x%x(C)" % (C1, row_major_index(0, 0, 4)*4+row_major_index(i, j, options.N)*16*4+offset_block_dense))
    print("  movaps %s, 0x%x(C)" % (C2, row_major_index(1, 0, 4)*4+row_major_index(i, j, options.N)*16*4+offset_block_dense))
    print("  movaps %s, 0x%x(C)" % (C3, row_major_index(2, 0, 4)*4+row_major_index(i, j, options.N)*16*4+offset_block_dense))
    print("  movaps %s, 0x%x(C)" % (C4, row_major_index(3, 0, 4)*4+row_major_index(i, j, options.N)*16*4+offset_block_dense))

    C1.release()
    C2.release()
    C3.release()
    C4.release()

# The main program.
parser = optparse.OptionParser(description =
"""This script generates a stream element kernel operating on 4x4 matrix
blocks. The kernel generated is written using assembly instructions assuming a
processor with the SSE instruction set.""")

parser.add_option("-N",
    metavar = "N",
    help = "generate kernel for NxN matrix of 4x4 matrix blocks [default: %default]",
    dest = "N",
    type = "int",
    default = 1)

parser.add_option("--name",
    metavar = "func",
    help = "set function name to \"func\" [default: %default]",
    dest = "functionName",
    type = "string",
    default = "stream_kernel")

parser.add_option("--no-checks",
    help = "generate code without any norm checks [default: %default]",
    action = "store_false",
    default = True,
    dest = "generate_checks")

parser.add_option("--no-alpha",
    help = "generate a kernel for the special case of alpha = 1 [default: %default]",
    action = "store_true",
    default = False,
    dest = "alphaOne")

parser.add_option("--debug",
    help = "print out a lot of debugging information [default: %default]",
    action = "store_true",
    default = False)

parser.add_option("--SSE",
    help = "set the SSE level from [1, 3, 4.1] [default: %default]",
    metavar = "level",
    type = "float",
    default = 1)

parser.add_option("--Z-curve",
    action = "store_true",
    default = False,
    help = """layout the multiply along a Z-curve as opposed to regular
row-major ordering [default: %default]""",
    dest = "Z_curve_ordering")

( options, arguments ) = parser.parse_args()

# Get logger.
log = logging.getLogger("generate assembly")
log.setLevel(logging.DEBUG)

logHandler = logging.StreamHandler()
if options.debug:
  logHandler.setLevel(logging.DEBUG)
else:
  logHandler.setLevel(logging.INFO)

logFormatter = logging.Formatter("%(levelname)s: %(message)s")
logHandler.setFormatter(logFormatter)
log.addHandler(logHandler)

# Check N.
if options.N <= 0:
  log.error("N needs to be a positive number > 0")
  sys.exit(1)

d = int(math.log(options.N)/math.log(2))
if 2**d != options.N:
  log.error("N needs to be a power of 2")
  sys.exit(1)

if options.Z_curve_ordering:
  if options.N != 4:
    log.warning("Z-curve ordering only works with N = 4")
    options.N = 4

# Check SSE level.
if not options.SSE in [1, 3, 4.1]:
  log.error("unknown SSE level")
  sys.exit(1)

# Generate assembly code.
print("# This code was generated by %s." % (os.path.basename(sys.argv[0])))
print("#")
print("# The command line was:")
print("#")
sys.stdout.write("# %s" % (os.path.basename(sys.argv[0])))
for i in range(1, len(sys.argv)):
  sys.stdout.write(" %s" % (sys.argv[i]))
sys.stdout.write("\n")
print("#")
print("# Code for a kernel matrix of %dx%d basic matrix blocks." % (options.N, options.N))

# Print some C function declarations.
print("")
print("# C API (as defined in spamm_kernel.h):")
print("#")
print("# struct spamm_multiply_stream_t")
print("# {")
print("#   struct spamm_data_t *A;")
print("#   struct spamm_data_t *B;")
print("#   struct spamm_data_t *C;")
print("# };")
print("#")
print("# void")
print("# %s (const unsigned int number_stream_elements," % (options.functionName))
print("#     float alpha,")
print("#     float tolerance,")
print("#     struct multiply_stream_t *multiply_stream);")
print("#")
print("# End of C API.")

print("")
print("# The matrix elements in the kernel block are layed out in the following order.")
print("# The basic matrix blocks are of size 4x4 to be able to take advantage fully of")
print("# single precision SSE instructions. Within the 4x4 blocks, the matrix elements")
print("# are layed out in row-major order. The blocks are themselves ordered in")
print("# row-major order within the kernel block.")

# Define some things.
print("")
print("# Function ABI.")
print("#define number_stream_elements %rdi")

# Get alpha if needed.
if not options.alphaOne:
  alpha = SSERegister("alpha", "%xmm0")
  print("#define %s %s" % (alpha.name, alpha.register))

# Get tolerance.
tolerance = SSERegister("tolerance", "%xmm1")

print("#define %s %s" % (tolerance.name, tolerance.register))
print("#define multiply_stream %rsi")

print("")
print("# Define loop variables.")
print("#define index %rax")
print("#define base_pointer %rdx")

print("")
print("# Define pointers to matrix data nodes in stream.")
print("#define A %r8")
print("#define B %r9")
print("#define C %r10")

if options.Z_curve_ordering:
  print("")
  print("# Define jump table variables.")
  print("#define jump_index %r11")
  print("#define jump_index_temp %r13")
  print("#define jump_index_base %rcx")
  print("#define jump_index_base_32 %ecx")
  print("#define old_stack %r12")

# The following sizes were generated with print_data_sizes.c.
sizeof_multiply_stream_t = 3*8

if not options.Z_curve_ordering:
  offset_norm = 16
  offset_block_dense = 192
  offset_block_dense_dilated = 1216
  offset_block_dense_transpose = 5312
else:
  offset_norm = 16
  offset_norm_upper = 192
  offset_norm_upper_transpose = 256
  offset_block_dense = 320
  offset_block_dense_dilated = 1344
  offset_block_dense_transpose = 5440

# Start the function prolog.
print("")
print("  # Function prolog.")
print("  .text")
print("  .balign 256")
print("  .global %s" % (options.functionName))
print("  .type %s, @function" % (options.functionName))

if options.Z_curve_ordering:
  print("")
  print("  # The jump table for the second kernel tier is in the read-only data section.")
  print("  .section .rodata")
  print("#ifdef __PIC__")
  print("  # For PIC we follow the example of gcc and do some fancier relative address")
  print("  # calculation to figure out the jump address for the jump table. Note that we")
  print("  # only need a 32-bit address (.long) although we compile this for 64-bit,")
  print("  # because of the relative address magic.")
  print("")
  print("  .balign 4")
  print("jump_table:")
  for jump_index in range(256):
    print("  .long tier_%03d-jump_table" % (jump_index))
  print("#else")
  print("  .balign 4")
  print("jump_table:")
  for jump_index in range(256):
    print("  .long tier_%03d" % (jump_index))
  print("#endif")
  print("  .text")

print("")
print("%s:" % (options.functionName))
print("")
print("  # Push used registers on stack.")
print("  push index")
print("  push base_pointer")
if options.Z_curve_ordering:
  print("  push jump_index")
  print("  push jump_index_temp")
  print("  push jump_index_base")
  print("  push old_stack")
print("  push A")
print("  push B")
print("  push C")

if options.Z_curve_ordering:
  print("")
  print("  # Push stack pointer so we can make room for local storage.")
  print("  mov %rsp, old_stack")
  print("  sub $0x16, %rsp")
  print("  and $-0x10, %rsp")
  print("")
  print("  # Load shuffle mask for norm comparisons.")
  print("  movl $0x0c080400, 0x0(%rsp)")
  print("  movl $0x80808080, 0x4(%rsp)")
  print("  movl $0x80808080, 0x8(%rsp)")
  print("  movl $0x80808080, 0xc(%rsp)")

  norm_mask = SSERegister("norm_mask")
  print("  movaps (%%rsp), %s" % (norm_mask))

if not options.alphaOne:
  print("")
  print("  # Copy alpha into all 4 elements of SSE register.")
  print("  shufps $0x0, alpha, alpha")

print("")
print("  # Test whether number_stream_elements is zero.")
print("  test number_stream_elements, number_stream_elements")
print("  jbe stream_done")
print("")
print("  # Set loop index to zero.")
print("  xor base_pointer, base_pointer")
print("  xor index, index")

block_counter = counter(1)

# Beginning of loop over stream elements.
print("")
print("  .balign 16")
print("stream_loop:")
print("")
print("  # Set the base pointer using sizeof(multiply_stream_t) = %d (0x%x)." % (sizeof_multiply_stream_t, sizeof_multiply_stream_t))
print("  imul $0x%x, base_pointer, base_pointer" % (sizeof_multiply_stream_t))
print("")
print("  # Load pointers to stream matrix blocks.")
print("  mov (multiply_stream, base_pointer, 1), A")
print("  mov 0x8(multiply_stream, base_pointer, 1), B")
print("  mov 0x10(multiply_stream, base_pointer, 1), C")

# Initialize the C registers so we can use them globally.
C1 = None
C2 = None
C3 = None
C4 = None

if options.Z_curve_ordering:
  print("")
  print("  # First level of hierarchy. Do some norm products [ A11*B11, A12*B21, A11*B12, A12*B22 ].")
  norm_1 = SSERegister("norm_1")
  print("  movaps 0x%x(A), %s" % (offset_norm_upper, norm_1))
  print("  mulps 0x%x(B), %s" % (offset_norm_upper_transpose, norm_1))
  print("  cmpps $0x02, %s, %s # norm product <= tolerance?" % (tolerance, norm_1))
  print("  pshufb %s, %s" % (norm_mask, norm_1))
  print("  pmovmskb %s, jump_index" % (norm_1))
  print("  shl $4, jump_index")
  print("")
  norm_2 = SSERegister("norm_2")
  print("  movaps 0x%x(A), %s" % (offset_norm_upper+4*4, norm_2))
  print("  mulps 0x%x(B), %s" % (offset_norm_upper_transpose+4*4, norm_2))
  print("  cmpps $0x02, %s, %s # norm product <= tolerance?" % (tolerance, norm_2))
  print("  pshufb %s, %s" % (norm_mask, norm_2))
  print("  pmovmskb %s, jump_index_temp" % (norm_2))
  print("")
  print("  or jump_index_temp, jump_index")
  norm_mask.release()
  norm_1.release()
  norm_2.release()

  print("")
  print("  # The value in jump_index is now such that a \"1\" bit indicates that the norm")
  print("  # product was <= tolerance and a \"0\" that it was not. The bits are ordered such")
  print("  # that:")
  print("  #")
  print("  # jump_index[0] = A11*B11")
  print("  # jump_index[1] = A12*B21")
  print("  # jump_index[2] = A11*B12")
  print("  # jump_index[3] = A12*B22")

  print("")
  print("  # Jump table for first tier.")
  print("#ifdef __PIC__")
  print("  # For PIC (Position Independent Code) we need to perform some magic")
  print("  # when it comes to figuring out the relative addresses of the jump")
  print("  # table. Fortunatley on x86_64 we can use %rdi.")
  print("")
  print("  lea jump_table(%rip), jump_index_base")
  print("  mov (jump_index_base, jump_index, 4), jump_index_base_32")
  print("  movslq jump_index_base_32, jump_index_base")
  print("  lea jump_table(%rip), jump_index")
  print("  lea (jump_index, jump_index_base), jump_index")
  print("  jmp *jump_index")
  print("#else")
  print("  jmp *jump_table(, jump_index, 4)")
  print("#endif")

  for jump_index in range(256):
    print("")
    print(".balign 16")
    print("tier_%03d:" % (jump_index))

    if (jump_index & (0x1 | 0x2)) == 0:
      print("  # Perform A_{11}*B_{11} = C_{11} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{11}*B_{11} + A_{12}*B_{21} = C_{11}")
      print("  # A_{11}*B_{12} + A_{12}*B_{22} = C_{12}")
      print("  # A_{21}*B_{11} + A_{22}*B_{21} = C_{21}")
      print("  # A_{21}*B_{12} + A_{22}*B_{22} = C_{22}")

      print("  # Perform A_{12}*B_{21} = C_{11} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{13}*B_{31} + A_{14}*B_{41} = C_{11}")
      print("  # A_{13}*B_{32} + A_{14}*B_{42} = C_{12}")
      print("  # A_{23}*B_{31} + A_{24}*B_{41} = C_{21}")
      print("  # A_{23}*B_{32} + A_{24}*B_{42} = C_{22}")

      block_product(1, 1, 1, clearC = True, writeC = False)
      block_product(1, 2, 1, clearC = False, writeC = False)
      block_product(1, 3, 1, clearC = False, writeC = False)
      block_product(1, 4, 1, clearC = False, writeC = True)
      block_product(1, 1, 2, clearC = True, writeC = False)
      block_product(1, 2, 2, clearC = False, writeC = False)
      block_product(1, 3, 2, clearC = False, writeC = False)
      block_product(1, 4, 2, clearC = False, writeC = True)
      block_product(2, 1, 1, clearC = True, writeC = False)
      block_product(2, 2, 1, clearC = False, writeC = False)
      block_product(2, 3, 1, clearC = False, writeC = False)
      block_product(2, 4, 1, clearC = False, writeC = True)
      block_product(2, 1, 2, clearC = True, writeC = False)
      block_product(2, 2, 2, clearC = False, writeC = False)
      block_product(2, 3, 2, clearC = False, writeC = False)
      block_product(2, 4, 2, clearC = False, writeC = True)

    else:
      if (jump_index & 0x1) == 0:
        print("  # Perform A_{11}*B_{11} = C_{11} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{11}*B_{11} + A_{12}*B_{21} = C_{11}")
        print("  # A_{11}*B_{12} + A_{12}*B_{22} = C_{12}")
        print("  # A_{21}*B_{11} + A_{22}*B_{21} = C_{21}")
        print("  # A_{21}*B_{12} + A_{22}*B_{22} = C_{22}")
        block_product(1, 1, 1, clearC = True, writeC = False)
        block_product(1, 2, 1, clearC = False, writeC = True)
        block_product(1, 1, 2, clearC = True, writeC = False)
        block_product(1, 2, 2, clearC = False, writeC = True)
        block_product(2, 1, 1, clearC = True, writeC = False)
        block_product(2, 2, 1, clearC = False, writeC = True)
        block_product(2, 1, 2, clearC = True, writeC = False)
        block_product(2, 2, 2, clearC = False, writeC = True)

      if (jump_index & 0x2) == 0:
        print("  # Perform A_{12}*B_{21} = C_{11} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{13}*B_{31} + A_{14}*B_{41} = C_{11}")
        print("  # A_{13}*B_{32} + A_{14}*B_{42} = C_{12}")
        print("  # A_{23}*B_{31} + A_{24}*B_{41} = C_{21}")
        print("  # A_{23}*B_{32} + A_{24}*B_{42} = C_{22}")
        block_product(1, 3, 1, clearC = True, writeC = False)
        block_product(1, 4, 1, clearC = False, writeC = True)
        block_product(1, 3, 2, clearC = True, writeC = False)
        block_product(1, 4, 2, clearC = False, writeC = True)
        block_product(2, 3, 1, clearC = True, writeC = False)
        block_product(2, 4, 1, clearC = False, writeC = True)
        block_product(2, 3, 2, clearC = True, writeC = False)
        block_product(2, 4, 2, clearC = False, writeC = True)

    if (jump_index & (0x4 | 0x8)) == 0:
      print("  # Perform A_{11}*B_{12} = C_{12} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{11}*B_{13} + A_{12}*B_{23} = C_{13}")
      print("  # A_{11}*B_{14} + A_{12}*B_{24} = C_{14}")
      print("  # A_{21}*B_{13} + A_{22}*B_{23} = C_{23}")
      print("  # A_{21}*B_{14} + A_{22}*B_{24} = C_{24}")

      print("  # Perform A_{12}*B_{22} = C_{12} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{13}*B_{33} + A_{14}*B_{43} = C_{13}")
      print("  # A_{13}*B_{34} + A_{14}*B_{44} = C_{14}")
      print("  # A_{23}*B_{33} + A_{24}*B_{43} = C_{23}")
      print("  # A_{23}*B_{34} + A_{24}*B_{44} = C_{24}")

      block_product(1, 1, 3, clearC = True, writeC = False)
      block_product(1, 2, 3, clearC = False, writeC = False)
      block_product(1, 3, 3, clearC = False, writeC = False)
      block_product(1, 4, 3, clearC = False, writeC = True)
      block_product(1, 1, 4, clearC = True, writeC = False)
      block_product(1, 2, 4, clearC = False, writeC = False)
      block_product(1, 3, 4, clearC = False, writeC = False)
      block_product(1, 4, 4, clearC = False, writeC = True)
      block_product(2, 1, 3, clearC = True, writeC = False)
      block_product(2, 2, 3, clearC = False, writeC = False)
      block_product(2, 3, 3, clearC = False, writeC = False)
      block_product(2, 4, 3, clearC = False, writeC = True)
      block_product(2, 1, 4, clearC = True, writeC = False)
      block_product(2, 2, 4, clearC = False, writeC = False)
      block_product(2, 3, 4, clearC = False, writeC = False)
      block_product(2, 4, 4, clearC = False, writeC = True)

    else:
      if (jump_index & 0x4) == 0:
        print("  # Perform A_{11}*B_{12} = C_{12} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{11}*B_{13} + A_{12}*B_{23} = C_{13}")
        print("  # A_{11}*B_{14} + A_{12}*B_{24} = C_{14}")
        print("  # A_{21}*B_{13} + A_{22}*B_{23} = C_{23}")
        print("  # A_{21}*B_{14} + A_{22}*B_{24} = C_{24}")
        block_product(1, 1, 3, clearC = True, writeC = False)
        block_product(1, 2, 3, clearC = False, writeC = True)
        block_product(1, 1, 4, clearC = True, writeC = False)
        block_product(1, 2, 4, clearC = False, writeC = True)
        block_product(2, 1, 3, clearC = True, writeC = False)
        block_product(2, 2, 3, clearC = False, writeC = True)
        block_product(2, 1, 4, clearC = True, writeC = False)
        block_product(2, 2, 4, clearC = False, writeC = True)

      if (jump_index & 0x8) == 0:
        print("  # Perform A_{12}*B_{22} = C_{12} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{13}*B_{33} + A_{14}*B_{43} = C_{13}")
        print("  # A_{13}*B_{34} + A_{14}*B_{44} = C_{14}")
        print("  # A_{23}*B_{33} + A_{24}*B_{43} = C_{23}")
        print("  # A_{23}*B_{34} + A_{24}*B_{44} = C_{24}")
        block_product(1, 3, 3, clearC = True, writeC = False)
        block_product(1, 4, 3, clearC = False, writeC = True)
        block_product(1, 3, 4, clearC = True, writeC = False)
        block_product(1, 4, 4, clearC = False, writeC = True)
        block_product(2, 3, 3, clearC = True, writeC = False)
        block_product(2, 4, 3, clearC = False, writeC = True)
        block_product(2, 3, 4, clearC = True, writeC = False)
        block_product(2, 4, 4, clearC = False, writeC = True)

    if (jump_index & (0x10 | 0x20)) == 0:
      print("  # Perform A_{21}*B_{11} = C_{21} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{31}*B_{11} + A_{32}*B_{21} = C_{31}")
      print("  # A_{31}*B_{12} + A_{32}*B_{22} = C_{32}")
      print("  # A_{41}*B_{11} + A_{42}*B_{21} = C_{41}")
      print("  # A_{41}*B_{12} + A_{42}*B_{22} = C_{42}")

      print("  # Perform A_{22}*B_{21} = C_{21} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{33}*B_{31} + A_{34}*B_{41} = C_{31}")
      print("  # A_{33}*B_{32} + A_{34}*B_{42} = C_{32}")
      print("  # A_{43}*B_{31} + A_{44}*B_{41} = C_{41}")
      print("  # A_{43}*B_{32} + A_{44}*B_{42} = C_{42}")

      block_product(3, 1, 1, clearC = True, writeC = False)
      block_product(3, 2, 1, clearC = False, writeC = False)
      block_product(3, 3, 1, clearC = False, writeC = False)
      block_product(3, 4, 1, clearC = False, writeC = True)
      block_product(3, 1, 2, clearC = True, writeC = False)
      block_product(3, 2, 2, clearC = False, writeC = False)
      block_product(3, 3, 2, clearC = False, writeC = False)
      block_product(3, 4, 2, clearC = False, writeC = True)
      block_product(4, 1, 1, clearC = True, writeC = False)
      block_product(4, 2, 1, clearC = False, writeC = False)
      block_product(4, 3, 1, clearC = False, writeC = False)
      block_product(4, 4, 1, clearC = False, writeC = True)
      block_product(4, 1, 2, clearC = True, writeC = False)
      block_product(4, 2, 2, clearC = False, writeC = False)
      block_product(4, 3, 2, clearC = False, writeC = False)
      block_product(4, 4, 2, clearC = False, writeC = True)

    else:
      if (jump_index & 0x10) == 0:
        print("  # Perform A_{21}*B_{11} = C_{21} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{31}*B_{11} + A_{32}*B_{21} = C_{31}")
        print("  # A_{31}*B_{12} + A_{32}*B_{22} = C_{32}")
        print("  # A_{41}*B_{11} + A_{42}*B_{21} = C_{41}")
        print("  # A_{41}*B_{12} + A_{42}*B_{22} = C_{42}")
        block_product(3, 1, 1, clearC = True, writeC = False)
        block_product(3, 2, 1, clearC = False, writeC = True)
        block_product(3, 1, 2, clearC = True, writeC = False)
        block_product(3, 2, 2, clearC = False, writeC = True)
        block_product(4, 1, 1, clearC = True, writeC = False)
        block_product(4, 2, 1, clearC = False, writeC = True)
        block_product(4, 1, 2, clearC = True, writeC = False)
        block_product(4, 2, 2, clearC = False, writeC = True)

      if (jump_index & 0x20) == 0:
        print("  # Perform A_{22}*B_{21} = C_{21} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{33}*B_{31} + A_{34}*B_{41} = C_{31}")
        print("  # A_{33}*B_{32} + A_{34}*B_{42} = C_{32}")
        print("  # A_{43}*B_{31} + A_{44}*B_{41} = C_{41}")
        print("  # A_{43}*B_{32} + A_{44}*B_{42} = C_{42}")
        block_product(3, 3, 1, clearC = True, writeC = False)
        block_product(3, 4, 1, clearC = False, writeC = True)
        block_product(3, 3, 2, clearC = True, writeC = False)
        block_product(3, 4, 2, clearC = False, writeC = True)
        block_product(4, 3, 1, clearC = True, writeC = False)
        block_product(4, 4, 1, clearC = False, writeC = True)
        block_product(4, 3, 2, clearC = True, writeC = False)
        block_product(4, 4, 2, clearC = False, writeC = True)

    if (jump_index & (0x30 | 0x40)) == 0:
      print("  # Perform A_{21}*B_{12} = C_{22} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{31}*B_{13} + A_{32}*B_{23} = C_{33}")
      print("  # A_{31}*B_{14} + A_{32}*B_{24} = C_{34}")
      print("  # A_{41}*B_{13} + A_{42}*B_{23} = C_{43}")
      print("  # A_{41}*B_{14} + A_{42}*B_{24} = C_{44}")

      print("  # Perform A_{22}*B_{22} = C_{22} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{33}*B_{33} + A_{34}*B_{43} = C_{33}")
      print("  # A_{33}*B_{34} + A_{34}*B_{44} = C_{34}")
      print("  # A_{43}*B_{33} + A_{44}*B_{43} = C_{43}")
      print("  # A_{43}*B_{34} + A_{44}*B_{44} = C_{44}")

      block_product(3, 1, 3, clearC = True, writeC = False)
      block_product(3, 2, 3, clearC = False, writeC = False)
      block_product(3, 3, 3, clearC = False, writeC = False)
      block_product(3, 4, 3, clearC = False, writeC = True)
      block_product(3, 1, 4, clearC = True, writeC = False)
      block_product(3, 2, 4, clearC = False, writeC = False)
      block_product(3, 3, 4, clearC = False, writeC = False)
      block_product(3, 4, 4, clearC = False, writeC = True)
      block_product(4, 1, 3, clearC = True, writeC = False)
      block_product(4, 2, 3, clearC = False, writeC = False)
      block_product(4, 3, 3, clearC = False, writeC = False)
      block_product(4, 4, 3, clearC = False, writeC = True)
      block_product(4, 1, 4, clearC = True, writeC = False)
      block_product(4, 2, 4, clearC = False, writeC = False)
      block_product(4, 3, 4, clearC = False, writeC = False)
      block_product(4, 4, 4, clearC = False, writeC = True)

    else:
      if (jump_index & 0x30) == 0:
        print("  # Perform A_{21}*B_{12} = C_{22} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{31}*B_{13} + A_{32}*B_{23} = C_{33}")
        print("  # A_{31}*B_{14} + A_{32}*B_{24} = C_{34}")
        print("  # A_{41}*B_{13} + A_{42}*B_{23} = C_{43}")
        print("  # A_{41}*B_{14} + A_{42}*B_{24} = C_{44}")
        block_product(3, 1, 3, clearC = True, writeC = False)
        block_product(3, 2, 3, clearC = False, writeC = True)
        block_product(3, 1, 4, clearC = True, writeC = False)
        block_product(3, 2, 4, clearC = False, writeC = True)
        block_product(4, 1, 3, clearC = True, writeC = False)
        block_product(4, 2, 3, clearC = False, writeC = True)
        block_product(4, 1, 4, clearC = True, writeC = False)
        block_product(4, 2, 4, clearC = False, writeC = True)

      if (jump_index & 0x40) == 0:
        print("  # Perform A_{22}*B_{22} = C_{22} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{33}*B_{33} + A_{34}*B_{43} = C_{33}")
        print("  # A_{33}*B_{34} + A_{34}*B_{44} = C_{34}")
        print("  # A_{43}*B_{33} + A_{44}*B_{43} = C_{43}")
        print("  # A_{43}*B_{34} + A_{44}*B_{44} = C_{44}")
        block_product(3, 3, 3, clearC = True, writeC = False)
        block_product(3, 4, 3, clearC = False, writeC = True)
        block_product(3, 3, 4, clearC = True, writeC = False)
        block_product(3, 4, 4, clearC = False, writeC = True)
        block_product(4, 3, 3, clearC = True, writeC = False)
        block_product(4, 4, 3, clearC = False, writeC = True)
        block_product(4, 3, 4, clearC = True, writeC = False)
        block_product(4, 4, 4, clearC = False, writeC = True)

else:
  for i in range(options.N):
    for j in range(options.N):
      for k in range(options.N):
        block_product(i+1, k+1, j+1, clearC = (k == 0), writeC = (k == options.N-1))

# End of outer loop.
print("")
print("loop_end:")
print("  # Loop end.")
print("  add $0x01, index")
print("  mov index, base_pointer")
print("  cmp  number_stream_elements, index")
print("  jb stream_loop")

# Leave function.
print("")
print("  .balign 16")
print("stream_done:")
print("")
print("  # Pop registers from stack.")
print("  pop C")
print("  pop B")
print("  pop A")
if options.Z_curve_ordering:
  print("  pop old_stack")
  print("  pop jump_index_base")
  print("  pop jump_index_temp")
  print("  pop jump_index")
print("  pop base_pointer")
print("  pop index")
print("")
print("  # Return from function.")
print("  ret")

if not options.alphaOne:
  alpha.release()
tolerance.release()

# Start function epilog.
print("")
print("  # Function epilog.")
print("  .size %s, .-%s" % (options.functionName, options.functionName))

# For debugging.
log.debug("remaining variables assigned as SSERegister: %s" % (SSERegister.variables))

if len(SSERegister.variables) > 0:
  log.error("still assigned variables as SSERegister: %s" % (SSERegister.variables))
  sys.exit(1)
