#!/usr/bin/python
#
# vim: tw=0

"""Generate SSE assembly code for a kernel operating on a 4x4 blocks.

The script generates a kernel written in assembly using SSE instructions that
operates on basic 4x4 matrix blocks. The kernel tier can be specified in the
function call.
"""

import math
import optparse
import os.path
import sys

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
processor with SSE2.""")

parser.add_option("-N",
    metavar = "N",
    help = "generate kernel for NxN matrix of 4x4 matrix blocks [default: %default]",
    dest = "N",
    type = "int",
    default = 1)

parser.add_option("--stripe",
    metavar = "N",
    help = "unroll loops over stripes consisting of N 4x4 blocks [default: %default]",
    dest = "N_stripe",
    type = "int",
    default = 1)

parser.add_option("--name",
    metavar = "func",
    help = "set function name to \"func\" [default: %default]",
    dest = "functionName",
    type = "string",
    default = "stream_kernel")

parser.add_option("--no-checks",
    action = "store_false",
    default = True,
    help = "generate code without any norm checks [default: %default]",
    dest = "generate_checks")

parser.add_option("--Z-curve",
    action = "store_true",
    default = False,
    help = """layout the multiply along a Z-curve as opposed to regular
row-major ordering [default: %default]""",
    dest = "Z_curve_ordering")

( options, arguments ) = parser.parse_args()

# Check N.
if options.N <= 0:
  print("N needs to be a positive number > 0")
  sys.exit(1)

d = int(math.log(options.N)/math.log(2))
if 2**d != options.N:
  print("N needs to be a power of 2")
  sys.exit(1)

# Check loop unrolling.
if options.N_stripe <= 0:
  print("unroll can not be zero or less, setting it to 1")
  options.N_stripe = 1

if options.N_stripe > options.N:
  print("unroll can not be greater than N, setting it to N (%d)" % (options.N))
  options.N_stripe = options.N

# Generate assembly code.
print("# This code was auto-generated by %s." % (os.path.basename(sys.argv[0])))
print("#")
print("# The command line was:")
print("#")
sys.stdout.write("# %s" % (os.path.basename(sys.argv[0])))
for i in range(1, len(sys.argv)):
  sys.stdout.write(" %s" % (sys.argv[i]))
sys.stdout.write("\n")
print("#")
print("# Code for a kernel matrix of %dx%d basic matrix blocks." % (options.N, options.N))
if options.N_stripe > 1:
  print("# Loops fully unrolled for a matrix of %dx%d basic matrix blocks." % (options.N_stripe, options.N_stripe))
else:
  print("# Loops not unrolled.")

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
print("# Define SSE registers used for C matrix")
print("#define C1 %xmm2")
print("#define C2 %xmm3")
print("#define C3 %xmm4")
print("#define C4 %xmm5")

print("")
print("# Define SSE registeres used for B matrix")
print("#define B1 %xmm6")
print("#define B2 %xmm7")
print("#define B3 %xmm8")
print("#define B4 %xmm9")

print("")
print("# Define SSE registeres used for A matrix")
print("#define A11 %xmm10")
print("#define A12 %xmm11")
print("#define A13 %xmm12")
print("#define A14 %xmm13")
print("#define A21 %xmm14")
print("#define A22 %xmm15")
print("#define A23 %xmm10")
print("#define A24 %xmm11")
print("#define A31 %xmm12")
print("#define A32 %xmm13")
print("#define A33 %xmm14")
print("#define A34 %xmm15")
print("#define A41 %xmm10")
print("#define A42 %xmm11")
print("#define A43 %xmm12")
print("#define A44 %xmm13")

print("")
print("# Define loop variables.")
print("#define index        %rax")
print("#define base_pointer %rdx")
if options.N-options.N_stripe > 0:
  print("#define i_index      %r11")
  print("#define j_index      %r12")

print("")
print("# Define pointers to matrix data nodes in stream.")
print("#define A %r8")
print("#define B %r9")
print("#define C %r10")

# The following sizes were generated with print_data_sizes.c.
sizeof_multiply_stream_t = 3*8
offset_norm = 24
offset_block_dense = 192
offset_block_dense_dilated = 1216

# Generate offsets.
print("")
print("# Define stream element size.")
print("#define SIZEOF_MULTIPLY_STREAM_T 0x%x" % (sizeof_multiply_stream_t))
print("")
print("# Define offsets into stream.")
print("#define OFFSET_NORM 0x%x" % (offset_norm))
print("#define OFFSET_BLOCK_DENSE 0x%x" % (offset_block_dense))
print("#define OFFSET_BLOCK_DENSE_DILATED 0x%x" % (offset_block_dense_dilated))
print("")
print("# Define offsets into A matrix.")
for i in range(options.N):
  for j in range(options.N):
    if options.Z_curve_ordering:
      print("#define A_OFFSET_%d%d %2d*64*4+OFFSET_BLOCK_DENSE_DILATED // %d = 0x%x" % (i+1, j+1,
          Z_curve_index(i, j),
          Z_curve_index(i, j)*64*4+offset_block_dense_dilated,
          Z_curve_index(i, j)*64*4+offset_block_dense_dilated))
    else:
      print("#define A_OFFSET_%d%d %2d*64*4+OFFSET_BLOCK_DENSE_DILATED // %d = 0x%x" % (i+1, j+1,
          row_major_index(i, j, options.N),
          row_major_index(i, j, options.N)*64*4+offset_block_dense_dilated,
          row_major_index(i, j, options.N)*64*4+offset_block_dense_dilated))

print("")
print("# Define offsets into B matrix.")
for i in range(options.N):
  for j in range(options.N):
    if options.Z_curve_ordering:
      print("#define B_OFFSET_%d%d %2d*16*4+OFFSET_BLOCK_DENSE // %d = 0x%x" % (i+1, j+1,
          Z_curve_index(i, j),
          Z_curve_index(i, j)*16*4+offset_block_dense,
          Z_curve_index(i, j)*16*4+offset_block_dense))
    else:
      print("#define B_OFFSET_%d%d %2d*16*4+OFFSET_BLOCK_DENSE // %d = 0x%x" % (i+1, j+1,
          row_major_index(i, j, options.N),
          row_major_index(i, j, options.N)*16*4+offset_block_dense,
          row_major_index(i, j, options.N)*16*4+offset_block_dense))

print("")
print("# Define offsets into C matrix.")
for i in range(options.N):
  for j in range(options.N):
    if options.Z_curve_ordering:
      print("#define C_OFFSET_%d%d %2d*16*4+OFFSET_BLOCK_DENSE // %d = 0x%x" % (i+1, j+1,
          Z_curve_index(i, j),
          Z_curve_index(i, j)*16*4+offset_block_dense,
          Z_curve_index(i, j)*16*4+offset_block_dense))
    else:
      print("#define C_OFFSET_%d%d %2d*16*4+OFFSET_BLOCK_DENSE // %d = 0x%x" % (i+1, j+1,
          row_major_index(i, j, options.N),
          row_major_index(i, j, options.N)*16*4+offset_block_dense,
          row_major_index(i, j, options.N)*16*4+offset_block_dense))

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
if options.N-options.N_stripe > 0:
  print("  push i_index")
  print("  push j_index")
print("  push A")
print("  push B")
print("  push C")

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
print("  .align 16")
print("stream_loop:")
print("")
print("  # Set the base pointer using sizeof(multiply_stream_t) = %d (0x%x)." % (sizeof_multiply_stream_t, sizeof_multiply_stream_t))
print("  imul $SIZEOF_MULTIPLY_STREAM_T, base_pointer, base_pointer")
print("")
print("  # Load pointers to stream matrix blocks.")
print("  mov (multiply_stream, base_pointer, 1), A")
print("  mov 0x8(multiply_stream, base_pointer, 1), B")
print("  mov 0x10(multiply_stream, base_pointer, 1), C")

# Loop over matrix blocks.
print("")
print("  # Loop over i index. Stream matrix block conists of %dx%d basic blocks." % (options.N, options.N))
print("  mov $%d, i_index" % (options.N))
print("")
print("  .align 16")
print("i_loop:")

print("")
print("  # Loop over j index. Stream matrix block conists of %dx%d basic blocks." % (options.N, options.N))
print("  mov $%d, j_index" % (options.N))
print("")
print("  .align 16")
print("j_loop:")

print("")
print("  # Adjust offset into matrix.")

for i in range(options.N_stripe):
  for j in range(options.N_stripe):
    print("")
    print("  # Reset C(%d,%d) matrix block accumulators." % (i+1, j+1))
    print("  xorps C1, C1")
    print("  xorps C2, C2")
    print("  xorps C3, C3")
    print("  xorps C4, C4")

    for k in range(options.N_stripe):
      if options.generate_checks:
        print("")
        print("  .align 16")
        print("jump_%d:" % (block_counter.get()))
        block_counter.increment()

        print("")
        print("  # Check norm of product ||A(%d,%d)||*||B(%d,%d)||." % (i+1, k+1, k+1, j+1))
        print("  movss 0x%x+OFFSET_NORM(A), B1" % ((i*options.N_stripe+k)*4))
        print("  mulss 0x%x+OFFSET_NORM(B), B1" % ((k*options.N_stripe+j)*4))
        print("  comiss tolerance, B1")
        print("  jb jump_%d" % (block_counter.get()))

      print("")
      print("  # Calculate C(%d,%d) += A(%d,%d)*B(%d,%d)." % (i+1, j+1, i+1, k+1, k+1, j+1))
      print("  movaps 0x%x+B_OFFSET_%d%d(B), B1" % (0*4*4, k+1, j+1))
      print("  movaps 0x%x+B_OFFSET_%d%d(B), B2" % (1*4*4, k+1, j+1))
      print("  movaps 0x%x+B_OFFSET_%d%d(B), B3" % (2*4*4, k+1, j+1))
      print("  movaps 0x%x+B_OFFSET_%d%d(B), B4" % (3*4*4, k+1, j+1))

      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(0, 0, 4)*4*4, i+1, k+1, 1, 1))
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(0, 1, 4)*4*4, i+1, k+1, 1, 2))
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(0, 2, 4)*4*4, i+1, k+1, 1, 3))
      print("  mulps B1, A11")
      print("  mulps B2, A12")
      print("  addps A11, C1")
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(0, 3, 4)*4*4, i+1, k+1, 1, 4))
      print("  mulps B3, A13")
      print("  addps A12, C1")
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(1, 0, 4)*4*4, i+1, k+1, 2, 1))
      print("  mulps B4, A14")
      print("  addps A13, C1")
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(1, 1, 4)*4*4, i+1, k+1, 2, 2))
      print("  mulps B1, A21")
      print("  addps A14, C1")
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(1, 2, 4)*4*4, i+1, k+1, 2, 3))
      print("  mulps B2, A22")
      print("  addps A21, C2")
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(1, 3, 4)*4*4, i+1, k+1, 2, 4))
      print("  mulps B3, A23")
      print("  addps A22, C2")
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(2, 0, 4)*4*4, i+1, k+1, 3, 1))
      print("  mulps B4, A24")
      print("  addps A23, C2")
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(2, 1, 4)*4*4, i+1, k+1, 3, 2))
      print("  mulps B1, A31")
      print("  addps A24, C2")
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(2, 2, 4)*4*4, i+1, k+1, 3, 3))
      print("  mulps B2, A32")
      print("  addps A31, C3")
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(2, 3, 4)*4*4, i+1, k+1, 3, 4))
      print("  mulps B3, A33")
      print("  addps A32, C3")
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(3, 0, 4)*4*4, i+1, k+1, 4, 1))
      print("  mulps B4, A34")
      print("  addps A33, C3")
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(3, 1, 4)*4*4, i+1, k+1, 4, 2))
      print("  mulps B1, A41")
      print("  addps A34, C3")
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(3, 2, 4)*4*4, i+1, k+1, 4, 3))
      print("  mulps B2, A42")
      print("  addps A41, C4")
      print("  movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(3, 3, 4)*4*4, i+1, k+1, 4, 4))
      print("  mulps B3, A43")
      print("  addps A42, C4")
      print("  mulps B4, A44")
      print("  addps A43, C4")
      print("  addps A44, C4")

    if options.generate_checks:
      print("")
      print("  .align 16")
      print("jump_%d:" % (block_counter.get()))
      block_counter.increment()

    print("")
    print("  # Multiply C(%d,%d) by alpha." % (i+1, j+1))
    print("  mulps alpha, C1")
    print("  mulps alpha, C2")
    print("  mulps alpha, C3")
    print("  mulps alpha, C4")

    print("")
    print("  # Add accumulated C(%d,%d) to already existing." % (i+1, j+1))
    print("  addps 0x0+C_OFFSET_%d%d(C), C1" % (i+1, j+1))
    print("  addps 0x10+C_OFFSET_%d%d(C), C2" % (i+1, j+1))
    print("  addps 0x20+C_OFFSET_%d%d(C), C3" % (i+1, j+1))
    print("  addps 0x30+C_OFFSET_%d%d(C), C4" % (i+1, j+1))

    print("")
    print("  # Write out C(%d,%d) submatrix block." % (i+1, j+1))
    print("  movaps C1, 0x0+C_OFFSET_%d%d(C)" % (i+1, j+1))
    print("  movaps C2, 0x10+C_OFFSET_%d%d(C)" % (i+1, j+1))
    print("  movaps C3, 0x20+C_OFFSET_%d%d(C)" % (i+1, j+1))
    print("  movaps C4, 0x30+C_OFFSET_%d%d(C)" % (i+1, j+1))

# End of inner loop.
if options.N-options.N_stripe > 0:
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
if options.N-options.N_stripe > 0:
  print("  pop j_index")
  print("  pop i_index")
print("  pop base_pointer")
print("  pop index")
print("")
print("  # Return from function.")
print("  ret")

# Start function epilog.
print("")
print("  # Function epilog.")
print("  .size %s, .-%s" % (options.functionName, options.functionName))
