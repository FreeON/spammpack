#!/usr/bin/python

"""Generate SSE assembly code for a kernel operating on a 4x4 blocks.

The script generates a kernel written in assembly using SSE instructions that
operates on basic 4x4 matrix blocks. The kernel tier can be specified in the
function call.
"""

import math, optparse, sys

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
    help = "generate fully unrolled kernel for NxN matrix of 4x4 matrix blocks [default: %default]",
    dest = "N",
    type = "int",
    default = 1)

parser.add_option("--unroll",
    metavar = "N",
    help = "fully unroll loops only at and below a matrix size of NxN [default: %default]",
    dest = "N_unroll",
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
  print "N needs to be a positive number > 0"
  sys.exit(1)

d = int(math.log(options.N)/math.log(2))
if 2**d != options.N:
  print "N needs to be a power of 2"
  sys.exit(1)

# Check loop unrolling.
if options.N_unroll <= 0:
  options.N_unroll = 1

if options.N_unroll > options.N:
  options.N_unroll = options.N

# Assembly code indentation.
padding = "  "

# Generate assembly code.
print "# This code was auto-generated by %s." % (sys.argv[0])
print "# The command line given was:"
print "#"
sys.stdout.write("# ")
for i in range(len(sys.argv)):
  sys.stdout.write(" %s" % (sys.argv[i]))
sys.stdout.write("\n")

# Print some C function declarations.
print "#"
print "# C API (found in spamm_kernel.h):"
print "#"
print "# struct spamm_multiply_stream_t"
print "# {"
print "#   struct spamm_data_t *A;"
print "#   struct spamm_data_t *B;"
print "#   struct spamm_data_t *C;"
print "# };"
print "#"
print "# void"
print "# %s (const unsigned int number_stream_elements," % (options.functionName)
print "#     float alpha,"
print "#     float tolerance,"
print "#     struct multiply_stream_t *multiply_stream);"
print "#"
print "# End of C API."

# Define some things.
print
print "# Function ABI."
print "#define number_stream_elements %rdi"
print "#define alpha                  %xmm0"
print "#define tolerance              %xmm1"
print "#define multiply_stream        %rsi"

print
print "# Define SSE registers used for C matrix"
print "#define C1 %xmm2"
print "#define C2 %xmm3"
print "#define C3 %xmm4"
print "#define C4 %xmm5"

print
print "# Define SSE registeres used for B matrix"
print "#define B1 %xmm6"
print "#define B2 %xmm7"
print "#define B3 %xmm8"
print "#define B4 %xmm9"

print
print "# Define SSE registeres used for A matrix"
print "#define A11 %xmm10"
print "#define A12 %xmm11"
print "#define A13 %xmm12"
print "#define A14 %xmm13"
print "#define A21 %xmm14"
print "#define A22 %xmm15"
print "#define A23 %xmm10"
print "#define A24 %xmm11"
print "#define A31 %xmm12"
print "#define A32 %xmm13"
print "#define A33 %xmm14"
print "#define A34 %xmm15"
print "#define A41 %xmm10"
print "#define A42 %xmm11"
print "#define A43 %xmm12"
print "#define A44 %xmm13"

print
print "# Define loop variables."
print "#define index        %rax"
print "#define base_pointer %rdx"

print
print "# Define pointers to matrix data nodes in stream."
print "#define A %r8"
print "#define B %r9"
print "#define C %r10"

sizeof_multiply_stream_t = 3*8
offset_norm = 24
offset_block_dense = 192
offset_block_dense_dilated = 1216

# Generate offsets.
print
print "# Define offsets into matrix blocks."
print
print "#define SIZEOF_MULTIPLY_STREAM_T 0x%x" % (sizeof_multiply_stream_t)
print
print "#define OFFSET_NORM 0x%x" % (offset_norm)
print "#define OFFSET_BLOCK_DENSE 0x%x" % (offset_block_dense)
print "#define OFFSET_BLOCK_DENSE_DILATED 0x%x" % (offset_block_dense_dilated)
print
for i in range(options.N):
  for j in range(options.N):
    if options.Z_curve_ordering:
      print "#define A_OFFSET_%d%d %2d*64*4+OFFSET_BLOCK_DENSE_DILATED // %d = 0x%x" % (i+1, j+1,
          Z_curve_index(i, j),
          Z_curve_index(i, j)*64*4+offset_block_dense_dilated,
          Z_curve_index(i, j)*64*4+offset_block_dense_dilated)
    else:
      print "#define A_OFFSET_%d%d %2d*64*4+OFFSET_BLOCK_DENSE_DILATED // %d = 0x%x" % (i+1, j+1,
          row_major_index(i, j, options.N),
          row_major_index(i, j, options.N)*64*4+offset_block_dense_dilated,
          row_major_index(i, j, options.N)*64*4+offset_block_dense_dilated)

print
for i in range(options.N):
  for j in range(options.N):
    if options.Z_curve_ordering:
      print "#define B_OFFSET_%d%d %2d*16*4+OFFSET_BLOCK_DENSE // %d = 0x%x" % (i+1, j+1,
          Z_curve_index(i, j),
          Z_curve_index(i, j)*16*4+offset_block_dense,
          Z_curve_index(i, j)*16*4+offset_block_dense)
    else:
      print "#define B_OFFSET_%d%d %2d*16*4+OFFSET_BLOCK_DENSE // %d = 0x%x" % (i+1, j+1,
          row_major_index(i, j, options.N),
          row_major_index(i, j, options.N)*16*4+offset_block_dense,
          row_major_index(i, j, options.N)*16*4+offset_block_dense)

print
for i in range(options.N):
  for j in range(options.N):
    if options.Z_curve_ordering:
      print "#define C_OFFSET_%d%d %2d*16*4+OFFSET_BLOCK_DENSE // %d = 0x%x" % (i+1, j+1,
          Z_curve_index(i, j),
          Z_curve_index(i, j)*16*4+offset_block_dense,
          Z_curve_index(i, j)*16*4+offset_block_dense)
    else:
      print "#define C_OFFSET_%d%d %2d*16*4+OFFSET_BLOCK_DENSE // %d = 0x%x" % (i+1, j+1,
          row_major_index(i, j, options.N),
          row_major_index(i, j, options.N)*16*4+offset_block_dense,
          row_major_index(i, j, options.N)*16*4+offset_block_dense)

# Start the function prolog.
print
print padding + "# Function prolog."
print padding + ".text"
print padding + ".align 256"
print padding + ".global %s" % (options.functionName)
print padding + ".type %s, @function" % (options.functionName)

print
print "%s:" % (options.functionName)
print
print padding + "# Push used registers on stack."
print padding + "push index"
print padding + "push base_pointer"
print padding + "push A"
print padding + "push B"
print padding + "push C"

print
print padding + "# Copy alpha into all 4 elements of SSE register."
print padding + "shufps $0x0, alpha, alpha"
print
print padding + "# Test whether number_stream_elements is zero."
print padding + "test number_stream_elements, number_stream_elements"
print padding + "jbe done"
print
print padding + "# Set loop index to zero."
print padding + "xor base_pointer, base_pointer"
print padding + "xor index, index"

block_counter = counter(1)

# Beginning of loop.
print
print padding + ".align 16"
print "loop:"
print
print padding + "# Set the base pointer using sizeof(multiply_stream_t) = %d (0x%x)." % (sizeof_multiply_stream_t, sizeof_multiply_stream_t)
print padding + "imul $SIZEOF_MULTIPLY_STREAM_T, base_pointer, base_pointer"
print
print padding + "# Load pointers to stream matrix blocks."
print padding + "mov (multiply_stream, base_pointer, 1), A"
print padding + "mov 0x8(multiply_stream, base_pointer, 1), B"
print padding + "mov 0x10(multiply_stream, base_pointer, 1), C"

for i in range(options.N):
  for j in range(options.N):
    print
    print padding + "# Reset C(%d,%d) matrix block accumulators." % (i+1, j+1)
    print padding + "xorps C1, C1"
    print padding + "xorps C2, C2"
    print padding + "xorps C3, C3"
    print padding + "xorps C4, C4"

    for k in range(options.N):
      if options.generate_checks:
        print
        print padding + ".align 16"
        print "block_%d:" % (block_counter.get())
        block_counter.increment()

        print
        print padding + "# Check norm of product ||A(%d,%d)||*||B(%d,%d)||." % (i+1, k+1, k+1, j+1)
        print padding + "movss 0x%x+OFFSET_NORM(A), B1" % ((i*options.N+k)*4)
        print padding + "mulss 0x%x+OFFSET_NORM(B), B1" % ((k*options.N+j)*4)
        print padding + "comiss tolerance, B1"
        print padding + "jb block_%d" % (block_counter.get())

      print
      print padding + "# Calculate C(%d,%d) += A(%d,%d)*B(%d,%d)." % (i+1, j+1, i+1, k+1, k+1, j+1)
      print padding + "movaps 0x%x+B_OFFSET_%d%d(B), B1" % (0*4*4, k+1, j+1)
      print padding + "movaps 0x%x+B_OFFSET_%d%d(B), B2" % (1*4*4, k+1, j+1)
      print padding + "movaps 0x%x+B_OFFSET_%d%d(B), B3" % (2*4*4, k+1, j+1)
      print padding + "movaps 0x%x+B_OFFSET_%d%d(B), B4" % (3*4*4, k+1, j+1)

      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(0, 0, 4)*4*4, i+1, k+1, 1, 1)
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(0, 1, 4)*4*4, i+1, k+1, 1, 2)
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(0, 2, 4)*4*4, i+1, k+1, 1, 3)
      print padding + "mulps B1, A11"
      print padding + "mulps B2, A12"
      print padding + "addps A11, C1"
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(0, 3, 4)*4*4, i+1, k+1, 1, 4)
      print padding + "mulps B3, A13"
      print padding + "addps A12, C1"
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(1, 0, 4)*4*4, i+1, k+1, 2, 1)
      print padding + "mulps B4, A14"
      print padding + "addps A13, C1"
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(1, 1, 4)*4*4, i+1, k+1, 2, 2)
      print padding + "mulps B1, A21"
      print padding + "addps A14, C1"
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(1, 2, 4)*4*4, i+1, k+1, 2, 3)
      print padding + "mulps B2, A22"
      print padding + "addps A21, C2"
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(1, 3, 4)*4*4, i+1, k+1, 2, 4)
      print padding + "mulps B3, A23"
      print padding + "addps A22, C2"
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(2, 0, 4)*4*4, i+1, k+1, 3, 1)
      print padding + "mulps B4, A24"
      print padding + "addps A23, C2"
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(2, 1, 4)*4*4, i+1, k+1, 3, 2)
      print padding + "mulps B1, A31"
      print padding + "addps A24, C2"
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(2, 2, 4)*4*4, i+1, k+1, 3, 3)
      print padding + "mulps B2, A32"
      print padding + "addps A31, C3"
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(2, 3, 4)*4*4, i+1, k+1, 3, 4)
      print padding + "mulps B3, A33"
      print padding + "addps A32, C3"
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(3, 0, 4)*4*4, i+1, k+1, 4, 1)
      print padding + "mulps B4, A34"
      print padding + "addps A33, C3"
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(3, 1, 4)*4*4, i+1, k+1, 4, 2)
      print padding + "mulps B1, A41"
      print padding + "addps A34, C3"
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(3, 2, 4)*4*4, i+1, k+1, 4, 3)
      print padding + "mulps B2, A42"
      print padding + "addps A41, C4"
      print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % (row_major_index(3, 3, 4)*4*4, i+1, k+1, 4, 4)
      print padding + "mulps B3, A43"
      print padding + "addps A42, C4"
      print padding + "mulps B4, A44"
      print padding + "addps A43, C4"
      print padding + "addps A44, C4"

    if options.generate_checks:
      print
      print padding + ".align 16"
      print "block_%d:" % (block_counter.get())
      block_counter.increment()

    print
    print padding + "# Multiply C(%d,%d) by alpha." % (i+1, j+1)
    print padding + "mulps alpha, C1"
    print padding + "mulps alpha, C2"
    print padding + "mulps alpha, C3"
    print padding + "mulps alpha, C4"

    print
    print padding + "# Add accumulated C(%d,%d) to already existing." % (i+1, j+1)
    print padding + "addps 0x0+C_OFFSET_%d%d(C), C1" % (i+1, j+1)
    print padding + "addps 0x10+C_OFFSET_%d%d(C), C2" % (i+1, j+1)
    print padding + "addps 0x20+C_OFFSET_%d%d(C), C3" % (i+1, j+1)
    print padding + "addps 0x30+C_OFFSET_%d%d(C), C4" % (i+1, j+1)

    print
    print padding + "# Write out C(%d,%d) submatrix block." % (i+1, j+1)
    print padding + "movaps C1, 0x0+C_OFFSET_%d%d(C)" % (i+1, j+1)
    print padding + "movaps C2, 0x10+C_OFFSET_%d%d(C)" % (i+1, j+1)
    print padding + "movaps C3, 0x20+C_OFFSET_%d%d(C)" % (i+1, j+1)
    print padding + "movaps C4, 0x30+C_OFFSET_%d%d(C)" % (i+1, j+1)

# End of loop.
print
print padding + "# Loop end."
print padding + "inc index"
print padding + "mov index, base_pointer"
print padding + "cmp number_stream_elements, index"
print padding + "jb loop"

# Leave function.
print
print padding + ".align 16"
print "done:"
print
print padding + "# Pop registers from stack."
print padding + "pop C"
print padding + "pop B"
print padding + "pop A"
print padding + "pop base_pointer"
print padding + "pop index"
print
print padding + "# Return from function."
print padding + "ret"

# Start function epilog.
print
print padding + "# Function epilog."
print padding + ".size %s, .-%s" % (options.functionName, options.functionName)
