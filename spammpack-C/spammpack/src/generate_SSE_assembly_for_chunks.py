#!/usr/bin/env python
#
# vim: tw=0

import argparse
import logging
import os.path
import sys
from SSERegister import SSERegister

#########################################################
#
# Global Constants.
#
#########################################################

sizeof_unsigned_int = 4

#########################################################
#
# Class Definitions.
#
#########################################################

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

#########################################################
#
# Function Definitions.
#
#########################################################

def row_major_index (i, j, N):
  """Return the index within a matrix block."""

  # Row-major storage.
  return i*N+j

def offset (i, j, N):
  """Return the offset into a matrix block."""

  return row_major_index(i, j, N)

def clearC (i, j):
  """Clears the C(i,j) accumulator registers."""

  # Fix up indices.
  i -= 1
  j -= 1

  # Reference the correct C registers.
  global C1
  global C2
  global C3
  global C4

  C1 = SSERegister(log, "C1")
  C2 = SSERegister(log, "C2")
  C3 = SSERegister(log, "C3")
  C4 = SSERegister(log, "C4")

  print("")
  print("  # Reset C(%d,%d) matrix block accumulators." % (i+1, j+1))
  print("  xorps %s, %s" % (C1, C1))
  print("  xorps %s, %s" % (C2, C2))
  print("  xorps %s, %s" % (C3, C3))
  print("  xorps %s, %s" % (C4, C4))

def writeC (i, j):
  """Write out the accumulator registers for C into memory."""

  # Fix up indices.
  i -= 1
  j -= 1

  # Reference the correct C registers.
  global C1
  global C2
  global C3
  global C4

  print("")
  print("  .balign 16")
  print("jump_%d:" % (block_counter.get()))
  block_counter.increment()

  print("")
  print("  # Multiply C(%d,%d) by alpha." % (i+1, j+1))
  print("  mulps alpha, %s" % (C1))
  print("  mulps alpha, %s" % (C2))
  print("  mulps alpha, %s" % (C3))
  print("  mulps alpha, %s" % (C4))

  print("")
  print("  # Add accumulated C(%d,%d) to already existing." % (i+1, j+1))
  print("  addps 0x%x(%%r15), %s" % (row_major_index(0, 0, 4)*4+offset(i, j, 4)*16*4, C1))
  print("  addps 0x%x(%%r15), %s" % (row_major_index(1, 0, 4)*4+offset(i, j, 4)*16*4, C2))
  print("  addps 0x%x(%%r15), %s" % (row_major_index(2, 0, 4)*4+offset(i, j, 4)*16*4, C3))
  print("  addps 0x%x(%%r15), %s" % (row_major_index(3, 0, 4)*4+offset(i, j, 4)*16*4, C4))

  print("")
  print("  # Write out C(%d,%d) submatrix block." % (i+1, j+1))
  print("  movaps %s, 0x%x(%%r15)" % (C1, row_major_index(0, 0, 4)*4+offset(i, j, 4)*16*4))
  print("  movaps %s, 0x%x(%%r15)" % (C2, row_major_index(1, 0, 4)*4+offset(i, j, 4)*16*4))
  print("  movaps %s, 0x%x(%%r15)" % (C3, row_major_index(2, 0, 4)*4+offset(i, j, 4)*16*4))
  print("  movaps %s, 0x%x(%%r15)" % (C4, row_major_index(3, 0, 4)*4+offset(i, j, 4)*16*4))

  C1.release()
  C2.release()
  C3.release()
  C4.release()

def block_product (i, k, j):
  """Produce an assembly code block to multiply 2 4x4 matrices in SSE. The
  index arguments are 1-based."""

  # Fix up indices.
  i -= 1
  j -= 1
  k -= 1

  # Reference the correct C registers.
  global C1
  global C2
  global C3
  global C4

  print("")
  print("  .balign 16")
  print("jump_%d:" % (block_counter.get()))
  block_counter.increment()

  norm = SSERegister(log, "norm")

  print("  # Check norm of product ||A(%d,%d)||*||B(%d,%d)||." % (i+1, k+1, k+1, j+1))
  print("#if SPAMM_NORM_TYPE == float")
  print("  movss 0x%x(%%r10), %s" % (row_major_index(i, k, 4)*4, norm))
  print("  mulss 0x%x(%%r11), %s" % (row_major_index(k, j, 4)*4, norm))

  # When comparing with the Intel Software Developer's Manual, keep in
  # mind that Intel uses Intel syntax and this code is writting using
  # At&T syntax, which means that the order of operand 1 and 2 are the
  # opposite which includes the comparison instruction.
  print("  comiss %s, %s" % (tolerance, norm))
  print("#elif SPAMM_NORM_TYPE == double")
  print("  movsd 0x%x(%%r10), %s" % (row_major_index(i, k, 4)*8, norm))
  print("  mulsd 0x%x(%%r11), %s" % (row_major_index(k, j, 4)*8, norm))
  print("  comisd %s, %s" % (tolerance, norm))
  print("#else")
  print("#error \"Unknown SPAMM_NORM_TYPE\"")
  print("#endif")
  print("  jbe jump_%d" % (block_counter.get()))

  norm.release()

  B1 = SSERegister(log, "B1")
  B2 = SSERegister(log, "B2")
  B3 = SSERegister(log, "B3")
  B4 = SSERegister(log, "B4")

  print("")
  print("  # Calculate C(%d,%d) += A(%d,%d)*B(%d,%d)." % (i+1, j+1, i+1, k+1, k+1, j+1))
  print("  movaps 0x%x(%%r14), %s" % (row_major_index(0, 0, 4)*4+offset(k, j, 4)*16*4, B1))
  print("  movaps 0x%x(%%r14), %s" % (row_major_index(1, 0, 4)*4+offset(k, j, 4)*16*4, B2))
  print("  movaps 0x%x(%%r14), %s" % (row_major_index(2, 0, 4)*4+offset(k, j, 4)*16*4, B3))
  print("  movaps 0x%x(%%r14), %s" % (row_major_index(3, 0, 4)*4+offset(k, j, 4)*16*4, B4))

  A11 = SSERegister(log, "A11")
  A12 = SSERegister(log, "A12")
  A13 = SSERegister(log, "A13")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(0, 0, 4)*4*4+offset(i, k, 4)*64*4, A11))
  print("  movaps 0x%x(%%r13), %s" % (row_major_index(0, 1, 4)*4*4+offset(i, k, 4)*64*4, A12))
  print("  movaps 0x%x(%%r13), %s" % (row_major_index(0, 2, 4)*4*4+offset(i, k, 4)*64*4, A13))
  print("  mulps %s, %s" % (B1, A11))
  print("  mulps %s, %s" % (B2, A12))
  print("  addps %s, %s" % (A11, C1))

  A11.release()
  A14 = SSERegister(log, "A14")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(0, 3, 4)*4*4+offset(i, k, 4)*64*4, A14))
  print("  mulps %s, %s" % (B3, A13))
  print("  addps %s, %s" % (A12, C1))

  A12.release()
  A21 = SSERegister(log, "A21")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(1, 0, 4)*4*4+offset(i, k, 4)*64*4, A21))
  print("  mulps %s, %s" % (B4, A14))
  print("  addps %s, %s" % (A13, C1))

  A13.release()
  A22 = SSERegister(log, "A22")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(1, 1, 4)*4*4+offset(i, k, 4)*64*4, A22))
  print("  mulps %s, %s" % (B1, A21))
  print("  addps %s, %s" % (A14, C1))

  A14.release()
  A23 = SSERegister(log, "A23")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(1, 2, 4)*4*4+offset(i, k, 4)*64*4, A23))
  print("  mulps %s, %s" % (B2, A22))
  print("  addps %s, %s" % (A21, C2))

  A21.release()
  A24 = SSERegister(log, "A24")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(1, 3, 4)*4*4+offset(i, k, 4)*64*4, A24))
  print("  mulps %s, %s" % (B3, A23))
  print("  addps %s, %s" % (A22, C2))

  A22.release()
  A31 = SSERegister(log, "A31")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(2, 0, 4)*4*4+offset(i, k, 4)*64*4, A31))
  print("  mulps %s, %s" % (B4, A24))
  print("  addps %s, %s" % (A23, C2))

  A23.release()
  A32 = SSERegister(log, "A32")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(2, 1, 4)*4*4+offset(i, k, 4)*64*4, A32))
  print("  mulps %s, %s" % (B1, A31))
  print("  addps %s, %s" % (A24, C2))

  A24.release()
  A33 = SSERegister(log, "A33")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(2, 2, 4)*4*4+offset(i, k, 4)*64*4, A33))
  print("  mulps %s, %s" % (B2, A32))
  print("  addps %s, %s" % (A31, C3))

  A31.release()
  A34 = SSERegister(log, "A34")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(2, 3, 4)*4*4+offset(i, k, 4)*64*4, A34))
  print("  mulps %s, %s" % (B3, A33))
  print("  addps %s, %s" % (A32, C3))

  A32.release()
  A41 = SSERegister(log, "A41")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(3, 0, 4)*4*4+offset(i, k, 4)*64*4, A41))
  print("  mulps %s, %s" % (B4, A34))
  print("  addps %s, %s" % (A33, C3))

  A33.release()
  A42 = SSERegister(log, "A42")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(3, 1, 4)*4*4+offset(i, k, 4)*64*4, A42))
  print("  mulps %s, %s" % (B1, A41))
  print("  addps %s, %s" % (A34, C3))

  A34.release()
  A43 = SSERegister(log, "A43")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(3, 2, 4)*4*4+offset(i, k, 4)*64*4, A43))
  print("  mulps %s, %s" % (B2, A42))
  print("  addps %s, %s" % (A41, C4))

  A41.release()
  A44 = SSERegister(log, "A44")

  print("  movaps 0x%x(%%r13), %s" % (row_major_index(3, 3, 4)*4*4+offset(i, k, 4)*64*4, A44))
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

#########################################################
#
# Main Program.
#
#########################################################

parser = argparse.ArgumentParser(description =
"""This script generates a stream element kernel operating on 4x4 matrix
blocks. The kernel generated is written using assembly instructions assuming a
processor with the SSE instruction set.""")

parser.add_argument("--name",
    metavar = "func",
    help = "set function name to \"func\" [default: %(default)s]",
    dest = "functionName",
    default = "spamm_stream_kernel")

parser.add_argument("--debug",
    help = "print out a lot of debugging information [default: %(default)s]",
    action = "store_true",
    default = False)

options = parser.parse_args()

# Setup logger.
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
print("# Code for a kernel matrix of %dx%d basic matrix blocks." % (4, 4))

# Print some C function declarations.
print("")
print("# C API (as defined in spamm_kernel.h):")
print("#")
print("# void")
print("# %s (const unsigned int number_stream_elements," % (options.functionName))
print("#     float alpha,")
print("#     float tolerance,")
print("#     unsigned int *stream,")
print("#     void *chunk_A,")
print("#     void *chunk_B,")
print("#     void *chunk_C)")
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
print("#")
print("# Following http://en.wikipedia.org/wiki/X86_calling_conventions#System_V_AMD64_ABI")
print("#")
print("# number_stream_elements -> %rdi")
print("# alpha                  -> %xmm0")
print("# tolerance              -> %xmm1")
print("# multiply_stream        -> %rsi")
print("# chunk_A                -> %rdx")
print("# chunk_B                -> %rcx")
print("# chunk_C                -> %r8")

print("")
print("#include \"config.h\"")

print("")
print("# Define some variables.")
print("#define number_stream_elements %rdi")

# Get alpha if needed.
alpha = SSERegister(log, "alpha", "%xmm0")
print("#define %s %s" % (alpha.name, alpha.register))

# Get tolerance.
tolerance = SSERegister(log, "tolerance", "%xmm1")

print("#define %s %s" % (tolerance.name, tolerance.register))
print("#define multiply_stream %rsi")
print("#define chunk_A         %rdx")
print("#define chunk_B         %rcx")
print("#define chunk_C         %r8")

print("")
print("# Define some constants for clarity of source.")
print("#define SIZEOF_INT 4")
print("#define SIZEOF_FLOAT 4")

print("")
print("# Define memory locations for spilling register.")
print("#define matrix_A_spill -0x08(%rbp)")
print("#define matrix_B_spill -0x10(%rbp)")
print("#define matrix_C_spill -0x18(%rbp)")
print("#define norm_A_spill   -0x20(%rbp)")
print("#define norm_B_spill   -0x28(%rbp)")
print("#define norm_C_spill   -0x30(%rbp)")

# Start the function prolog.
print("")
print("  # Function prolog.")
print("  .text")
print("  .balign 256")
print("  .global %s" % (options.functionName))
print("  .type %s, @function" % (options.functionName))
print("")
print("%s:" % (options.functionName))

print("")
print("  # Save a few registers.")
print("  push %rax")
print("  push %rbx")
print("  push %r9")
print("  push %r10")
print("  push %r11")
print("  push %r12")
print("  push %r13")
print("  push %r14")
print("  push %r15")

print("")
print("  # Create some local storage.")
print("  push %rbp")
print("  mov %rsp, %rbp")
print("  sub $0x40, %rsp")

print("")
print("  # Copy alpha into all 4 elements of SSE register.")
print("  shufps $0x0, alpha, alpha")
print("")
print("  # Test whether number_stream_elements is zero.")
print("  test number_stream_elements, number_stream_elements")
print("  jbe stream_done")

block_counter = counter(1)

# Some preparations.
print("")
print("  # Load pointers to stream matrix blocks.")
print("  mov 6*8(chunk_A), %r10 # A_dilated.")
print("  mov 5*8(chunk_B), %r11 # B.")
print("  mov 5*8(chunk_C), %r12 # C.")
print("  lea (chunk_A, %r10), %r10")
print("  lea (chunk_B, %r11), %r11")
print("  lea (chunk_C, %r12), %r12")
print("  mov %r10, matrix_A_spill")
print("  mov %r11, matrix_B_spill")
print("  mov %r12, matrix_C_spill")

print("")
print("  # Calculate offset into norm arrays.")
print("  movl 1*4(chunk_A), %ebx  # Number of tiers.")
print("  sub $1, %rbx             # number_tiers-1.")
print("  xor %rax, %rax           # Loop counter.")
print("  xor %r9, %r9             # Offset into norm array.")
print("  mov $1, %r10             # Load the number of norms of tier 0.")
print("")
print("  .balign 16")
print("tier_loop:")
print("  add %r10, %r9")
print("  sal $2, %r10             # %r10 <- 4*%r10.")
print("  add $1, %rax")
print("  cmp %rbx, %rax")
print("  jne tier_loop")

print("")
print("  # Load pointers to tier norms.")
print("  mov 7*8(chunk_A), %r10")
print("  mov 7*8(chunk_B), %r11")
print("  mov 7*8(chunk_C), %r12")
print("  lea (chunk_A, %r10), %r10")
print("  lea (chunk_B, %r11), %r11")
print("  lea (chunk_C, %r12), %r12")
print("  lea (%r10, %r9, SIZEOF_FLOAT), %r10")
print("  lea (%r11, %r9, SIZEOF_FLOAT), %r11")
print("  lea (%r12, %r9, SIZEOF_FLOAT), %r12")
print("  mov %r10, norm_A_spill")
print("  mov %r11, norm_B_spill")
print("  mov %r12, norm_C_spill")

print("")
print("  # Beginning of loop over stream elements.")

print("")
print("  # Set loop index to zero.")
print("  xor %rax, %rax")

print("")
print("  .balign 16")
print("stream_loop:")

print("")
print("  # Set the base pointer to the correct offset in the multiply_stream.")
print("  imul $3, %rax, %r8")
print("")
print("  # Load the linear indices of the next 16x16 block. Note that since we use")
print("  # unsigned int for the linear index, a 32-bit value, we need to load a")
print("  # long (32-bit) value.")
print("  movl  (multiply_stream, %r8, SIZEOF_INT), %r10d")
print("  movl 4(multiply_stream, %r8, SIZEOF_INT), %r11d")
print("  movl 8(multiply_stream, %r8, SIZEOF_INT), %r12d")
print("  mov %r10, %r13")
print("  mov %r11, %r14")
print("  mov %r12, %r15")

print("")
print("  # Calculate base of norm offset. SPAMM_N_KERNEL x SPAMM_N_KERNEL blocks")
print("  # are Z-curve ordered, and there is one norm per SPAMM_N_BLOCK x SPAMM_N_BLOCK")
print("  # matrix.")
print("  imul $SPAMM_N_KERNEL_BLOCKED*SPAMM_N_KERNEL_BLOCKED*SIZEOF_FLOAT, %r10")
print("  imul $SPAMM_N_KERNEL_BLOCKED*SPAMM_N_KERNEL_BLOCKED*SIZEOF_FLOAT, %r11")
print("  imul $SPAMM_N_KERNEL_BLOCKED*SPAMM_N_KERNEL_BLOCKED*SIZEOF_FLOAT, %r12")
print("  add norm_A_spill, %r10")
print("  add norm_B_spill, %r11")
print("  add norm_C_spill, %r12")

print("")
print("  # Calculate base of matrix offset.")
print("  imul $SPAMM_N_KERNEL*SPAMM_N_KERNEL*SIZEOF_FLOAT*4, %r13  # Dilated matrix A.")
print("  imul $SPAMM_N_KERNEL*SPAMM_N_KERNEL*SIZEOF_FLOAT,   %r14  # Matrix B.")
print("  imul $SPAMM_N_KERNEL*SPAMM_N_KERNEL*SIZEOF_FLOAT,   %r15  # Matrix C.")
print("  add matrix_A_spill, %r13")
print("  add matrix_B_spill, %r14")
print("  add matrix_C_spill, %r15")

# Initialize the C registers so we can use them globally.
C1 = None
C2 = None
C3 = None
C4 = None

for i in range(4):
  for j in range(4):
    clearC(i+1, j+1)
    for k in range(4):
      block_product(i+1, k+1, j+1)
    writeC(i+1, j+1)

# End of outer loop.
print("")
print("loop_end:")
print("  # Loop end.")
print("  add $0x01, %rax")
print("  cmp  number_stream_elements, %rax")
print("  jb stream_loop")

# Leave function.
print("")
print("  .balign 16")
print("stream_done:")

print("")
print("  # Restore old stack frame.")
print("  leave")

print("")
print("  # Get registers back from stack.")
print("  pop %r15")
print("  pop %r14")
print("  pop %r13")
print("  pop %r12")
print("  pop %r11")
print("  pop %r10")
print("  pop %r9")
print("  pop %rbx")
print("  pop %rax")

print("")
print("  # Return from function.")
print("  ret")

alpha.release()
tolerance.release()

# Start function epilog.
print("")
print("  # Function epilog.")
print("  .size %s, .-%s" % (options.functionName, options.functionName))

if len(SSERegister.variables) > 0:
  print("still assigned variables as SSERegister: %s" % (SSERegister.variables))
  sys.exit(1)
