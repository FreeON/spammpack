#!/usr/bin/python
#
# Generate SSE assembly code for a kernel operating on a 4x4 blocks.

import math, optparse, sys

class box:
  def __init__ (self, i_1, i_2, j_1, j_2):
    self.i_1 = i_1
    self.i_2 = i_2
    self.j_1 = j_1
    self.j_2 = j_2

  def __str__ (self):
    return "box: [%d-%d][%d-%d]" % (self.i_1, self.i_2, self.j_1, self.j_2)

class counter:
  def __init__ (self):
    self.counter = 0

  def __init__ (self, initial_value):
    self.counter = initial_value

  def increment (self):
    self.counter += 1

  def get (self):
    return self.counter

# Generate matrix product with Z-curve ordering.
def generate_Z_curve (A, B, C, block_counter):
  if A.i_2-A.i_1 == 1 and A.j_2-A.j_1 == 1:
    i = C.i_1
    j = C.j_1
    k = A.j_1

    if options.generate_checks:
      print
      print padding + ".align 16"
      print "block_%d:" % (block_counter.get())
      block_counter.increment()

      print
      print padding + "# Check norm of product A(%d,%d)*B(%d,%d)." % (i+1, k+1, k+1, j+1)
      print padding + "movss 0x%x(multiply_stream, base_pointer), B1" % ((i*options.N+k)*4+24)
      print padding + "mulss 0x%x(multiply_stream, base_pointer), B1" % ((k*options.N+j+options.N**2)*4+24)
      print padding + "comiss tolerance, B1"
      print padding + "jb block_%d" % (block_counter.get())

    print
    print padding + "# Reset C(%d,%d) matrix block accumulators." % (i+1, j+1)
    print padding + "xorps C1, C1"
    print padding + "xorps C2, C2"
    print padding + "xorps C3, C3"
    print padding + "xorps C4, C4"

    print
    print padding + "# Calculate C(%d,%d) = A(%d,%d)*B(%d,%d)." % (i+1, j+1, i+1, k+1, k+1, j+1)
    print padding + "movaps 0x%x+B_OFFSET_%d%d(B), B1" % (0*4*4, k+1, j+1)
    print padding + "movaps 0x%x+B_OFFSET_%d%d(B), B2" % (1*4*4, k+1, j+1)
    print padding + "movaps 0x%x+B_OFFSET_%d%d(B), B3" % (2*4*4, k+1, j+1)
    print padding + "movaps 0x%x+B_OFFSET_%d%d(B), B4" % (3*4*4, k+1, j+1)

    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((0*4+0)*4*4, i+1, k+1, 1, 1)
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((0*4+1)*4*4, i+1, k+1, 1, 2)
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((0*4+2)*4*4, i+1, k+1, 1, 3)
    print padding + "mulps B1, A11"
    print padding + "mulps B2, A12"
    print padding + "addps A11, C1"
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((0*4+3)*4*4, i+1, k+1, 1, 4)
    print padding + "mulps B3, A13"
    print padding + "addps A12, C1"
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((1*4+0)*4*4, i+1, k+1, 2, 1)
    print padding + "mulps B4, A14"
    print padding + "addps A13, C1"
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((1*4+1)*4*4, i+1, k+1, 2, 2)
    print padding + "mulps B1, A21"
    print padding + "addps A14, C1"
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((1*4+2)*4*4, i+1, k+1, 2, 3)
    print padding + "mulps B2, A22"
    print padding + "addps A21, C2"
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((1*4+3)*4*4, i+1, k+1, 2, 4)
    print padding + "mulps B3, A23"
    print padding + "addps A22, C2"
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((2*4+0)*4*4, i+1, k+1, 3, 1)
    print padding + "mulps B4, A24"
    print padding + "addps A23, C2"
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((2*4+1)*4*4, i+1, k+1, 3, 2)
    print padding + "mulps B1, A31"
    print padding + "addps A24, C2"
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((2*4+2)*4*4, i+1, k+1, 3, 3)
    print padding + "mulps B2, A32"
    print padding + "addps A31, C3"
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((2*4+3)*4*4, i+1, k+1, 3, 4)
    print padding + "mulps B3, A33"
    print padding + "addps A32, C3"
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((3*4+0)*4*4, i+1, k+1, 4, 1)
    print padding + "mulps B4, A34"
    print padding + "addps A33, C3"
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((3*4+1)*4*4, i+1, k+1, 4, 2)
    print padding + "mulps B1, A41"
    print padding + "addps A34, C3"
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((3*4+2)*4*4, i+1, k+1, 4, 3)
    print padding + "mulps B2, A42"
    print padding + "addps A41, C4"
    print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((3*4+3)*4*4, i+1, k+1, 4, 4)
    print padding + "mulps B3, A43"
    print padding + "addps A42, C4"
    print padding + "mulps B4, A44"
    print padding + "addps A43, C4"
    print padding + "addps A44, C4"

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

  else:
    A_11 = box(A.i_1, A.i_1+(A.i_2-A.i_1)/2, A.j_1, A.j_1+(A.j_2-A.j_1)/2)
    A_12 = box(A.i_1, A.i_1+(A.i_2-A.i_1)/2, A.j_1+(A.j_2-A.j_1)/2, A.j_2)
    A_21 = box(A.i_1+(A.i_2-A.i_1)/2, A.i_2, A.j_1, A.j_1+(A.j_2-A.j_1)/2)
    A_22 = box(A.i_1+(A.i_2-A.i_1)/2, A.i_2, A.j_1+(A.j_2-A.j_1)/2, A.j_2)

    B_11 = box(B.i_1, B.i_1+(B.i_2-B.i_1)/2, B.j_1, B.j_1+(B.j_2-B.j_1)/2)
    B_12 = box(B.i_1, B.i_1+(B.i_2-B.i_1)/2, B.j_1+(B.j_2-B.j_1)/2, B.j_2)
    B_21 = box(B.i_1+(B.i_2-B.i_1)/2, B.i_2, B.j_1, B.j_1+(B.j_2-B.j_1)/2)
    B_22 = box(B.i_1+(B.i_2-B.i_1)/2, B.i_2, B.j_1+(B.j_2-B.j_1)/2, B.j_2)

    C_11 = box(C.i_1, C.i_1+(C.i_2-C.i_1)/2, C.j_1, C.j_1+(C.j_2-C.j_1)/2)
    C_12 = box(C.i_1, C.i_1+(C.i_2-C.i_1)/2, C.j_1+(C.j_2-C.j_1)/2, C.j_2)
    C_21 = box(C.i_1+(C.i_2-C.i_1)/2, C.i_2, C.j_1, C.j_1+(C.j_2-C.j_1)/2)
    C_22 = box(C.i_1+(C.i_2-C.i_1)/2, C.i_2, C.j_1+(C.j_2-C.j_1)/2, C.j_2)

    generate_Z_curve(A_11, B_11, C_11, block_counter)
    generate_Z_curve(A_12, B_21, C_11, block_counter)
    generate_Z_curve(A_11, B_12, C_12, block_counter)
    generate_Z_curve(A_12, B_22, C_12, block_counter)
    generate_Z_curve(A_21, B_11, C_21, block_counter)
    generate_Z_curve(A_22, B_21, C_21, block_counter)
    generate_Z_curve(A_21, B_12, C_22, block_counter)
    generate_Z_curve(A_22, B_22, C_22, block_counter)

# Main program.
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
#print "#define i_outer      %r10"
#print "#define j_outer      %r11"

print
print "# Define pointers to matrix blocks in stream."
print "#define A %r8"
print "#define B %rcx"
print "#define C %r9"

# Generate offsets.
print
print "# Define offsets into matrix blocks."
print
for i in range(options.N):
  for j in range(options.N):
    print "#define A_OFFSET_%d%d (%d*%d+%d)*64*4 // %d = 0x%x" % (i+1, j+1, i, options.N, j, (i*options.N+j)*64, (i*options.N+j)*64)

print
for i in range(options.N):
  for j in range(options.N):
    print "#define B_OFFSET_%d%d (%d*%d+%d)*16*4 // %d = 0x%x" % (i+1, j+1, i, options.N, j, (i*options.N+j)*16, (i*options.N+j)*16)

print
for i in range(options.N):
  for j in range(options.N):
    print "#define C_OFFSET_%d%d (%d*%d+%d)*16*4 // %d = 0x%x" % (i+1, j+1, i, options.N, j, (i*options.N+j)*16, (i*options.N+j)*16)

# Print some C function declarations.
print
print "# C function declaration"
print "#"
print "# struct multiply_stream_t"
print "# {"
print "#   float *A_block;"
print "#   float *B_block;"
print "#   float *C_block;"
print "#   float  norm[%d];" % (2*options.N**2)
print "# };"
print "#"
print "# void"
print "# %s (const unsigned int number_stream_elements," % (options.functionName)
print "#     float alpha,"
print "#     float tolerance,"
print "#     struct multiply_stream_t *multiply_stream);"

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
#print padding + "push i_outer"
#print padding + "push j_outer"

print
print padding + "# Copy alpha into all 4 elements of SSE register."
print padding + "shufps $0x0, alpha, alpha"
print
print padding + "# Divide number of stream elements by %d to simulate stride of %d." % (options.N**3, options.N**3)
print padding + "shr $%i, number_stream_elements" % (3*math.log(options.N)/math.log(2))
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
print padding + "# Set the base pointer using sizeof(multiply_stream_t) = 0x98."
print padding + "imul $0x98, base_pointer, base_pointer"
print
print padding + "# Load pointers to stream matrix blocks."
print padding + "mov (multiply_stream, base_pointer, 1), A"
print padding + "mov 0x8(multiply_stream, base_pointer, 1), B"
print padding + "mov 0x10(multiply_stream, base_pointer, 1), C"

if options.Z_curve_ordering:
  generate_Z_curve(box(0, options.N, 0, options.N),
      box(0, options.N, 0, options.N),
      box(0, options.N, 0, options.N),
      block_counter)

  if options.generate_checks:
    print
    print padding + ".align 16"
    print "block_%d:" % (block_counter.get())
    block_counter.increment()

else:
  #if options.N_unroll < options.N:
  #  # Generate outer loop code.
  #  print
  #  print padding + ".align 16"
  #  print "outer_i:"

  #if options.N_unroll < options.N:
  #  # Generate outer loop code.
  #  print
  #  print padding + ".align 16"
  #  print "outer_j:"

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
          print padding + "# Check norm of product."
          print padding + "movss 0x%x(multiply_stream, base_pointer), B1" % ((i*options.N+k)*4+24)
          print padding + "mulss 0x%x(multiply_stream, base_pointer), B1" % ((k*options.N+j+options.N**2)*4+24)
          print padding + "comiss tolerance, B1"
          print padding + "jb block_%d" % (block_counter.get())

        print
        print padding + "# Calculate C(%d,%d) = A(%d,%d)*B(%d,%d)." % (i+1, j+1, i+1, k+1, k+1, j+1)
        print padding + "movaps 0x%x+B_OFFSET_%d%d(B), B1" % (0*4*4, k+1, j+1)
        print padding + "movaps 0x%x+B_OFFSET_%d%d(B), B2" % (1*4*4, k+1, j+1)
        print padding + "movaps 0x%x+B_OFFSET_%d%d(B), B3" % (2*4*4, k+1, j+1)
        print padding + "movaps 0x%x+B_OFFSET_%d%d(B), B4" % (3*4*4, k+1, j+1)

        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((0*4+0)*4*4, i+1, k+1, 1, 1)
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((0*4+1)*4*4, i+1, k+1, 1, 2)
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((0*4+2)*4*4, i+1, k+1, 1, 3)
        print padding + "mulps B1, A11"
        print padding + "mulps B2, A12"
        print padding + "addps A11, C1"
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((0*4+3)*4*4, i+1, k+1, 1, 4)
        print padding + "mulps B3, A13"
        print padding + "addps A12, C1"
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((1*4+0)*4*4, i+1, k+1, 2, 1)
        print padding + "mulps B4, A14"
        print padding + "addps A13, C1"
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((1*4+1)*4*4, i+1, k+1, 2, 2)
        print padding + "mulps B1, A21"
        print padding + "addps A14, C1"
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((1*4+2)*4*4, i+1, k+1, 2, 3)
        print padding + "mulps B2, A22"
        print padding + "addps A21, C2"
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((1*4+3)*4*4, i+1, k+1, 2, 4)
        print padding + "mulps B3, A23"
        print padding + "addps A22, C2"
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((2*4+0)*4*4, i+1, k+1, 3, 1)
        print padding + "mulps B4, A24"
        print padding + "addps A23, C2"
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((2*4+1)*4*4, i+1, k+1, 3, 2)
        print padding + "mulps B1, A31"
        print padding + "addps A24, C2"
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((2*4+2)*4*4, i+1, k+1, 3, 3)
        print padding + "mulps B2, A32"
        print padding + "addps A31, C3"
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((2*4+3)*4*4, i+1, k+1, 3, 4)
        print padding + "mulps B3, A33"
        print padding + "addps A32, C3"
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((3*4+0)*4*4, i+1, k+1, 4, 1)
        print padding + "mulps B4, A34"
        print padding + "addps A33, C3"
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((3*4+1)*4*4, i+1, k+1, 4, 2)
        print padding + "mulps B1, A41"
        print padding + "addps A34, C3"
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((3*4+2)*4*4, i+1, k+1, 4, 3)
        print padding + "mulps B2, A42"
        print padding + "addps A41, C4"
        print padding + "movaps 0x%x+A_OFFSET_%d%d(A), A%d%d" % ((3*4+3)*4*4, i+1, k+1, 4, 4)
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
#print padding + "pop j_outer"
#print padding + "pop i_outer"
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
