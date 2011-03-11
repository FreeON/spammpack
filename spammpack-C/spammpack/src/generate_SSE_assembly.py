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

class register:
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

    if name in register.variables:
      log.error("The variable %s is already assigned a register" % (name))
      sys.exit(1)

    if registerName and not registerName in register.registerPool:
      log.error("The requested register %s is not available anymore" % (registerName))
      sys.exit(1)

    if registerName:
      self.name = name
      self.register = registerName
      register.registerPool.remove(registerName)
      register.variables[name] = registerName
    else:
      if len(register.registerPool) > 0:
        self.name = name
        self.register = register.registerPool.pop(0)
        register.variables[name] = self.register
      else:
        log.error("no registers left")
        sys.exit(1)

    log.debug("assigned %s --> %s" % (self.register, self.name))
    log.debug("registerPool: %s" % (register.registerPool))
    log.debug("variables: %s" % (register.variables))

  def release (self):
    """Release a register back to the register pool."""

    if not self.name in register.variables:
      log.error("The variable %s is not assigned a register" % (self.name))
      sys.exit(1)

    register.registerPool.append(self.register)
    del register.variables[self.name]

    log.debug("released register %s assigned to variable %s" % (self.register, self.name))
    log.debug("registerPool: %s" % (register.registerPool))
    log.debug("variables: %s" % (register.variables))

  def __str__ (self):
    """Use this variable in the code."""

    if not self.name in register.variables:
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
    help = "set the SSE level from [1, 4.1] [default: %default]",
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

# Check SSE level.
if not options.SSE in [1, 4.1]:
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
print("#define alpha                  %xmm0")
print("#define tolerance              %xmm1")
print("#define multiply_stream        %rsi")

print("")
print("# Define loop variables.")
print("#define index        %rax")
print("#define base_pointer %rdx")

print("")
print("# Define pointers to matrix data nodes in stream.")
print("#define A %r8")
print("#define B %r9")
print("#define C %r10")

# The following sizes were generated with print_data_sizes.c.
sizeof_multiply_stream_t = 3*8
offset_norm = 16
offset_block_dense = 192
offset_block_dense_dilated = 1216
offset_block_dense_transpose = 5312

# Start the function prolog.
print("")
print("  # Function prolog.")
print("  .text")
print("  .align 256")
print("  .global %s" % (options.functionName))
print("  .type %s, @function" % (options.functionName))

print("")
print("%s:" % (options.functionName))
print("")
print("  # Push used registers on stack.")
print("  push index")
print("  push base_pointer")
print("  push A")
print("  push B")
print("  push C")

if not options.alphaOne:
  print("")
  alpha = register("alpha", "%xmm0")
  print("  # Copy alpha into all 4 elements of SSE register.")
  print("  shufps $0x0, alpha, alpha")

# Get tolerance.
tolerance = register("tolerance", "%xmm1")

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
print("  .align 16")
print("stream_loop:")
print("")
print("  # Set the base pointer using sizeof(multiply_stream_t) = %d (0x%x)." % (sizeof_multiply_stream_t, sizeof_multiply_stream_t))
print("  imul $0x%x, base_pointer, base_pointer" % (sizeof_multiply_stream_t))
print("")
print("  # Load pointers to stream matrix blocks.")
print("  mov (multiply_stream, base_pointer, 1), A")
print("  mov 0x8(multiply_stream, base_pointer, 1), B")
print("  mov 0x10(multiply_stream, base_pointer, 1), C")

for i in range(options.N):
  for j in range(options.N):

    C1 = register("C1")
    C2 = register("C2")
    C3 = register("C3")
    C4 = register("C4")

    print("")
    print("  # Reset C(%d,%d) matrix block accumulators." % (i+1, j+1))
    print("  xorps %s, %s" % (C1, C1))
    print("  xorps %s, %s" % (C2, C2))
    print("  xorps %s, %s" % (C3, C3))
    print("  xorps %s, %s" % (C4, C4))

    for k in range(options.N):
      if options.generate_checks:
        print("")
        print("  .align 16")
        print("jump_%d:" % (block_counter.get()))
        block_counter.increment()

        norm = register("norm")

        print("")
        print("  # Check norm of product ||A(%d,%d)||*||B(%d,%d)||." % (i+1, k+1, k+1, j+1))
        print("  movss 0x%x(A), %s" % ((i*options.N+k)*4+offset_norm, norm))
        print("  mulss 0x%x(B), %s" % ((k*options.N+j)*4+offset_norm, norm))

        # When comparing with the Intel Software Developer's Manual, keep in
        # mind that Intel uses Intel syntax and this code is writting using
        # At&T syntax, which means that the order of operand 1 and 2 are the
        # opposite.
        print("  comiss tolerance, %s" % (norm))
        print("  jbe jump_%d" % (block_counter.get()))

        norm.release()

      if options.SSE == 1:
        B1 = register("B1")
        B2 = register("B2")
        B3 = register("B3")
        B4 = register("B4")

        print("")
        print("  # Calculate C(%d,%d) += A(%d,%d)*B(%d,%d)." % (i+1, j+1, i+1, k+1, k+1, j+1))
        print("  movaps 0x%x(B), %s" % (row_major_index(0, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense, B1))
        print("  movaps 0x%x(B), %s" % (row_major_index(1, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense, B2))
        print("  movaps 0x%x(B), %s" % (row_major_index(2, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense, B3))
        print("  movaps 0x%x(B), %s" % (row_major_index(3, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense, B4))

        A11 = register("A11")
        A12 = register("A12")
        A13 = register("A13")

        print("  movaps 0x%x(A), %s" % (row_major_index(0, 0, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A11))
        print("  movaps 0x%x(A), %s" % (row_major_index(0, 1, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A12))
        print("  movaps 0x%x(A), %s" % (row_major_index(0, 2, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A13))
        print("  mulps %s, %s" % (B1, A11))
        print("  mulps %s, %s" % (B2, A12))
        print("  addps %s, %s" % (A11, C1))

        A11.release()
        A14 = register("A14")

        print("  movaps 0x%x(A), %s" % (row_major_index(0, 3, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A14))
        print("  mulps %s, %s" % (B3, A13))
        print("  addps %s, %s" % (A12, C1))

        A12.release()
        A21 = register("A21")

        print("  movaps 0x%x(A), %s" % (row_major_index(1, 0, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A21))
        print("  mulps %s, %s" % (B4, A14))
        print("  addps %s, %s" % (A13, C1))

        A13.release()
        A22 = register("A22")

        print("  movaps 0x%x(A), %s" % (row_major_index(1, 1, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A22))
        print("  mulps %s, %s" % (B1, A21))
        print("  addps %s, %s" % (A14, C1))

        A14.release()
        A23 = register("A23")

        print("  movaps 0x%x(A), %s" % (row_major_index(1, 2, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A23))
        print("  mulps %s, %s" % (B2, A22))
        print("  addps %s, %s" % (A21, C2))

        A21.release()
        A24 = register("A24")

        print("  movaps 0x%x(A), %s" % (row_major_index(1, 3, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A24))
        print("  mulps %s, %s" % (B3, A23))
        print("  addps %s, %s" % (A22, C2))

        A22.release()
        A31 = register("A31")

        print("  movaps 0x%x(A), %s" % (row_major_index(2, 0, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A31))
        print("  mulps %s, %s" % (B4, A24))
        print("  addps %s, %s" % (A23, C2))

        A23.release()
        A32 = register("A32")

        print("  movaps 0x%x(A), %s" % (row_major_index(2, 1, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A32))
        print("  mulps %s, %s" % (B1, A31))
        print("  addps %s, %s" % (A24, C2))

        A24.release()
        A33 = register("A33")

        print("  movaps 0x%x(A), %s" % (row_major_index(2, 2, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A33))
        print("  mulps %s, %s" % (B2, A32))
        print("  addps %s, %s" % (A31, C3))

        A31.release()
        A34 = register("A34")

        print("  movaps 0x%x(A), %s" % (row_major_index(2, 3, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A34))
        print("  mulps %s, %s" % (B3, A33))
        print("  addps %s, %s" % (A32, C3))

        A32.release()
        A41 = register("A41")

        print("  movaps 0x%x(A), %s" % (row_major_index(3, 0, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A41))
        print("  mulps %s, %s" % (B4, A34))
        print("  addps %s, %s" % (A33, C3))

        A33.release()
        A42 = register("A42")

        print("  movaps 0x%x(A), %s" % (row_major_index(3, 1, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A42))
        print("  mulps %s, %s" % (B1, A41))
        print("  addps %s, %s" % (A34, C3))

        A34.release()
        A43 = register("A43")

        print("  movaps 0x%x(A), %s" % (row_major_index(3, 2, 4)*4*4+row_major_index(i, k, options.N)*64*4+offset_block_dense_dilated, A43))
        print("  mulps %s, %s" % (B2, A42))
        print("  addps %s, %s" % (A41, C4))

        A41.release()
        A44 = register("A44")

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

      elif options.SSE == 4.1:
        B1 = register("B1")
        B2 = register("B2")
        B3 = register("B3")
        B4 = register("B4")

        print("")
        print("  # Calculate C(%d,%d) += A(%d,%d)*B(%d,%d)." % (i+1, j+1, i+1, k+1, k+1, j+1))
        print("  movaps 0x%x(B), %s" % (row_major_index(0, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense_transpose, B1))
        print("  movaps 0x%x(B), %s" % (row_major_index(1, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense_transpose, B2))
        print("  movaps 0x%x(B), %s" % (row_major_index(2, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense_transpose, B3))
        print("  movaps 0x%x(B), %s" % (row_major_index(3, 0, 4)*4+row_major_index(k, j, options.N)*16*4+offset_block_dense_transpose, B4))

        A1 = register("A1")
        A2 = register("A2")

        print("  movaps 0x%x(A), %s" % (row_major_index(0, 0, 4)*4+row_major_index(i, k, options.N)*16*4+offset_block_dense, A1))
        print("  movaps 0x%x(A), %s" % (row_major_index(1, 0, 4)*4+row_major_index(i, k, options.N)*16*4+offset_block_dense, A2))

        print("")
        print("  # Calculate C(1,:).")
        C11 = register("C11")

        print("  movaps %s, %s" % (B1, C11))
        print("  dpps $0xf1, %s, %s" % (A1, C11))

        C12 = register("C12")

        print("  movaps %s, %s" % (B2, C12))
        print("  dpps $0xf2, %s, %s" % (A1, C12))

        C13 = register("C13")

        print("  movaps %s, %s" % (B3, C13))
        print("  dpps $0xf4, %s, %s" % (A1, C13))

        C14 = register("C14")

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

        A3 = register("A3")
        print("  movaps 0x%x(A), %s" % (row_major_index(2, 0, 4)*4+row_major_index(i, k, options.N)*16*4+offset_block_dense, A3))

        print("")
        print("  # Calculate C(2,:).")
        C21 = register("C21")

        print("  movaps %s, %s" % (B1, C21))
        print("  dpps $0xf1, %s, %s" % (A2, C21))

        C22 = register("C22")

        print("  movaps %s, %s" % (B2, C22))
        print("  dpps $0xf2, %s, %s" % (A2, C22))

        C23 = register("C23")

        print("  movaps %s, %s" % (B3, C23))
        print("  dpps $0xf4, %s, %s" % (A2, C23))

        C24 = register("C24")

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

        A4 = register("A4")
        print("  movaps 0x%x(A), %s" % (row_major_index(3, 0, 4)*4+row_major_index(i, k, options.N)*16*4+offset_block_dense, A4))

        print("")
        print("  # Calculate C(3,:).")
        C31 = register("C31")

        print("  movaps %s, %s" % (B1, C31))
        print("  dpps $0xf1, %s, %s" % (A3, C31))

        C32 = register("C32")

        print("  movaps %s, %s" % (B2, C32))
        print("  dpps $0xf2, %s, %s" % (A3, C32))

        C33 = register("C33")

        print("  movaps %s, %s" % (B3, C33))
        print("  dpps $0xf4, %s, %s" % (A3, C33))

        C34 = register("C34")

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

        print("")
        print("  # Calculate C(4,:).")
        C41 = register("C41")

        print("  movaps %s, %s" % (B1, C41))
        print("  dpps $0xf1, %s, %s" % (A4, C41))

        C42 = register("C42")

        print("  movaps %s, %s" % (B2, C42))
        print("  dpps $0xf2, %s, %s" % (A4, C42))

        C43 = register("C43")

        print("  movaps %s, %s" % (B3, C43))
        print("  dpps $0xf4, %s, %s" % (A4, C43))

        C44 = register("C44")

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

    if options.generate_checks:
      print("")
      print("  .align 16")
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

# End of inner loop.
if options.N-options.N > 0:
  print("")
  print("  # End of loop over j index.")
  print("  sub $%d, A" % (0))
  print("  sub $%d, B" % (0))
  print("  sub $%d, C" % (0))
  print("  dec j_index")
  print("  test j_index, j_index")
  print("  ja j_loop")

  print("")
  print("  # End of loop over i index.")
  print("  dec i_index")
  print("  test i_index, i_index")
  print("  ja i_loop")

# End of outer loop.
print("")
print("  # Loop end.")
print("  inc index")
print("  mov index, base_pointer")
print("  cmp number_stream_elements, index")
print("  jb stream_loop")

# Leave function.
print("")
print("  .align 16")
print("stream_done:")
print("")
print("  # Pop registers from stack.")
print("  pop C")
print("  pop B")
print("  pop A")
if options.N-options.N > 0:
  print("  pop j_index")
  print("  pop i_index")
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
