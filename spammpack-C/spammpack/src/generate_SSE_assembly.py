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
import argparse
import os.path
import sys
from spammOffsets import spammOffsets
from SSERegister import SSERegister

class Line:
  """A class that helps print a line of assembly code. It takes care of things
  such as line numbers and intendation."""

  def __init__ (self):
    if Line.initialized:
      log.error("Only one Line class can be used at a time")
      sys.exit(1)

    Line.initialized = True
    self.line_number = 0

  def out (self, line):
    """Print a line."""

    self.line_number += 1
    print(line)

  def get_line_number (self):
    """Get the current line number (the line number of the last line
    printed)."""

    return self.line_number

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

def offset (i, j, N):
  """Return the offset into a matrix block."""

  if options.Z_curve_ordering:
    return Z_curve_index(i, j)
  else:
    return row_major_index(i, j, N)

def check_load_offset (load_offset):
  """Check the overlap with the load offset and the last store offset."""

  global last_store_offset

  for store_offset in last_store_offset:
    log.debug("checking load offset 0x%x against store offset 0x%x, overlap is 0x%x"
        % (load_offset, store_offset, abs(store_offset-load_offset)%4096))
    if abs(store_offset-load_offset)%4096 == 0 and abs(store_offset-load_offset) > 0:
      log.warn("4 kB overlap detected, load offset 0x%x against store offset 0x%x, diff = 0x%x, overlap is 0x%x"
          % (load_offset, store_offset, abs(store_offset-load_offset), abs(store_offset-load_offset)%4096))

def issue_load (load_offset, base_pointer, register):
  """Issue a load statement from memory to a register. The memory address is
  given by load_offset(base_pointer)."""

  check_load_offset(load_offset)
  print("  movaps 0x%x(%s), %s" % (load_offset, base_pointer, register))

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

def block_product (i, k, j, number_deactivated_products):
  """Produce an assembly code block to multiply 2 4x4 matrices in SSE. The
  index arguments are 1-based."""

  # Check whether to omit this product.
  if number_deactivated_products > 0:
    number_deactivated_products -= 1
    return number_deactivated_products

  # Fix up indices.
  i -= 1
  j -= 1
  k -= 1

  # Reference the correct C registers.
  global C1
  global C2
  global C3
  global C4

  global last_store_offset

  if options.generate_checks:
    print("")
    print("  .balign 16")
    print("jump_%d:" % (block_counter.get()))
    block_counter.increment()

    norm = SSERegister(log, "norm")

    print("  # Check norm of product ||A(%d,%d)||*||B(%d,%d)||." % (i+1, k+1, k+1, j+1))
    check_load_offset((i*options.N+k)*4+spammOffsets.offset_norm)
    print("  movss 0x%x(A), %s" % ((i*options.N+k)*4+spammOffsets.offset_norm, norm))
    check_load_offset((k*options.N+j)*4+spammOffsets.offset_norm)
    print("  mulss 0x%x(B), %s" % ((k*options.N+j)*4+spammOffsets.offset_norm, norm))

    # When comparing with the Intel Software Developer's Manual, keep in
    # mind that Intel uses Intel syntax and this code is writting using
    # At&T syntax, which means that the order of operand 1 and 2 are the
    # opposite which includes the comparison instruction.
    print("  comiss %s, %s" % (tolerance, norm))
    print("  jbe jump_%d" % (block_counter.get()))

    norm.release()

  if options.SSE == 1:
    B1 = SSERegister(log, "B1")
    B2 = SSERegister(log, "B2")
    B3 = SSERegister(log, "B3")
    B4 = SSERegister(log, "B4")

    print("")
    print("  # Calculate C(%d,%d) += A(%d,%d)*B(%d,%d)." % (i+1, j+1, i+1, k+1, k+1, j+1))
    print("  movaps 0x%x(B), %s" % (row_major_index(0, 0, 4)*4+offset(k, j, options.N)*16*4+spammOffsets.offset_block_dense, B1))
    print("  movaps 0x%x(B), %s" % (row_major_index(1, 0, 4)*4+offset(k, j, options.N)*16*4+spammOffsets.offset_block_dense, B2))
    print("  movaps 0x%x(B), %s" % (row_major_index(2, 0, 4)*4+offset(k, j, options.N)*16*4+spammOffsets.offset_block_dense, B3))
    print("  movaps 0x%x(B), %s" % (row_major_index(3, 0, 4)*4+offset(k, j, options.N)*16*4+spammOffsets.offset_block_dense, B4))

    A11 = SSERegister(log, "A11")
    A12 = SSERegister(log, "A12")
    A13 = SSERegister(log, "A13")

    print("  movaps 0x%x(A), %s" % (row_major_index(0, 0, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A11))
    print("  movaps 0x%x(A), %s" % (row_major_index(0, 1, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A12))
    print("  movaps 0x%x(A), %s" % (row_major_index(0, 2, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A13))
    print("  mulps %s, %s" % (B1, A11))
    print("  mulps %s, %s" % (B2, A12))
    print("  addps %s, %s" % (A11, C1))

    A11.release()
    A14 = SSERegister(log, "A14")

    print("  movaps 0x%x(A), %s" % (row_major_index(0, 3, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A14))
    print("  mulps %s, %s" % (B3, A13))
    print("  addps %s, %s" % (A12, C1))

    A12.release()
    A21 = SSERegister(log, "A21")

    print("  movaps 0x%x(A), %s" % (row_major_index(1, 0, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A21))
    print("  mulps %s, %s" % (B4, A14))
    print("  addps %s, %s" % (A13, C1))

    A13.release()
    A22 = SSERegister(log, "A22")

    print("  movaps 0x%x(A), %s" % (row_major_index(1, 1, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A22))
    print("  mulps %s, %s" % (B1, A21))
    print("  addps %s, %s" % (A14, C1))

    A14.release()
    A23 = SSERegister(log, "A23")

    print("  movaps 0x%x(A), %s" % (row_major_index(1, 2, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A23))
    print("  mulps %s, %s" % (B2, A22))
    print("  addps %s, %s" % (A21, C2))

    A21.release()
    A24 = SSERegister(log, "A24")

    print("  movaps 0x%x(A), %s" % (row_major_index(1, 3, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A24))
    print("  mulps %s, %s" % (B3, A23))
    print("  addps %s, %s" % (A22, C2))

    A22.release()
    A31 = SSERegister(log, "A31")

    print("  movaps 0x%x(A), %s" % (row_major_index(2, 0, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A31))
    print("  mulps %s, %s" % (B4, A24))
    print("  addps %s, %s" % (A23, C2))

    A23.release()
    A32 = SSERegister(log, "A32")

    print("  movaps 0x%x(A), %s" % (row_major_index(2, 1, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A32))
    print("  mulps %s, %s" % (B1, A31))
    print("  addps %s, %s" % (A24, C2))

    A24.release()
    A33 = SSERegister(log, "A33")

    print("  movaps 0x%x(A), %s" % (row_major_index(2, 2, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A33))
    print("  mulps %s, %s" % (B2, A32))
    print("  addps %s, %s" % (A31, C3))

    A31.release()
    A34 = SSERegister(log, "A34")

    print("  movaps 0x%x(A), %s" % (row_major_index(2, 3, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A34))
    print("  mulps %s, %s" % (B3, A33))
    print("  addps %s, %s" % (A32, C3))

    A32.release()
    A41 = SSERegister(log, "A41")

    print("  movaps 0x%x(A), %s" % (row_major_index(3, 0, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A41))
    print("  mulps %s, %s" % (B4, A34))
    print("  addps %s, %s" % (A33, C3))

    A33.release()
    A42 = SSERegister(log, "A42")

    print("  movaps 0x%x(A), %s" % (row_major_index(3, 1, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A42))
    print("  mulps %s, %s" % (B1, A41))
    print("  addps %s, %s" % (A34, C3))

    A34.release()
    A43 = SSERegister(log, "A43")

    print("  movaps 0x%x(A), %s" % (row_major_index(3, 2, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A43))
    print("  mulps %s, %s" % (B2, A42))
    print("  addps %s, %s" % (A41, C4))

    A41.release()
    A44 = SSERegister(log, "A44")

    print("  movaps 0x%x(A), %s" % (row_major_index(3, 3, 4)*4*4+offset(i, k, options.N)*64*4+spammOffsets.offset_block_dense_dilated, A44))
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

    A1 = SSERegister(log, "A1")
    issue_load(row_major_index(0, 0, 4)*4+offset(i, k, options.N)*16*4+spammOffsets.offset_block_dense, "A", A1)

    B1 = SSERegister(log, "B1")
    issue_load(row_major_index(0, 0, 4)*4+offset(k, j, options.N)*16*4+spammOffsets.offset_block_dense_transpose, "B", B1)
    B2 = SSERegister(log, "B2")
    issue_load(row_major_index(1, 0, 4)*4+offset(k, j, options.N)*16*4+spammOffsets.offset_block_dense_transpose, "B", B2)
    B3 = SSERegister(log, "B3")
    issue_load(row_major_index(2, 0, 4)*4+offset(k, j, options.N)*16*4+spammOffsets.offset_block_dense_transpose, "B", B3)
    B4 = SSERegister(log, "B4")
    issue_load(row_major_index(3, 0, 4)*4+offset(k, j, options.N)*16*4+spammOffsets.offset_block_dense_transpose, "B", B4)

    print("")
    print("  # Calculate C(1,:).")

    C11 = SSERegister(log, "C11")
    print("  movaps %s, %s" % (B1, C11))
    print("  dpps $0xf1, %s, %s" % (A1, C11))

    C12 = SSERegister(log, "C12")
    print("  movaps %s, %s" % (B2, C12))
    print("  dpps $0xf2, %s, %s" % (A1, C12))
    print("  blendps $0x01, %s, %s" % (C11, C12))
    C11.release()

    C13 = SSERegister(log, "C13")
    print("  movaps %s, %s" % (B3, C13))
    print("  dpps $0xf4, %s, %s" % (A1, C13))
    print("  blendps $0x03, %s, %s" % (C12, C13))
    C12.release()

    C14 = SSERegister(log, "C14")
    print("  movaps %s, %s" % (B4, C14))
    print("  dpps $0xf8, %s, %s" % (A1, C14))
    print("  blendps $0x07, %s, %s" % (C13, C14))
    C13.release()

    print("  addps %s, %s" % (C14, C1))
    C14.release()
    A1.release()

    A2 = SSERegister(log, "A2")
    issue_load(row_major_index(1, 0, 4)*4+offset(i, k, options.N)*16*4+spammOffsets.offset_block_dense, "A", A2)

    print("")
    print("  # Calculate C(2,:).")

    C21 = SSERegister(log, "C21")
    print("  movaps %s, %s" % (B1, C21))
    print("  dpps $0xf1, %s, %s" % (A2, C21))

    C22 = SSERegister(log, "C22")
    print("  movaps %s, %s" % (B2, C22))
    print("  dpps $0xf2, %s, %s" % (A2, C22))
    print("  blendps $0x01, %s, %s" % (C21, C22))
    C21.release()

    C23 = SSERegister(log, "C23")
    print("  movaps %s, %s" % (B3, C23))
    print("  dpps $0xf4, %s, %s" % (A2, C23))
    print("  blendps $0x03, %s, %s" % (C22, C23))
    C22.release()

    C24 = SSERegister(log, "C24")
    print("  movaps %s, %s" % (B4, C24))
    print("  dpps $0xf8, %s, %s" % (A2, C24))
    print("  blendps $0x07, %s, %s" % (C23, C24))
    C23.release()

    print("  addps %s, %s" % (C24, C2))
    C24.release()
    A2.release()

    A3 = SSERegister(log, "A3")
    issue_load(row_major_index(2, 0, 4)*4+offset(i, k, options.N)*16*4+spammOffsets.offset_block_dense, "A", A3)

    print("")
    print("  # Calculate C(3,:).")

    C31 = SSERegister(log, "C31")
    print("  movaps %s, %s" % (B1, C31))
    print("  dpps $0xf1, %s, %s" % (A3, C31))

    C32 = SSERegister(log, "C32")
    print("  movaps %s, %s" % (B2, C32))
    print("  dpps $0xf2, %s, %s" % (A3, C32))
    print("  blendps $0x01, %s, %s" % (C31, C32))
    C31.release()

    C33 = SSERegister(log, "C33")
    print("  movaps %s, %s" % (B3, C33))
    print("  dpps $0xf4, %s, %s" % (A3, C33))
    print("  blendps $0x03, %s, %s" % (C32, C33))
    C32.release()

    C34 = SSERegister(log, "C34")
    print("  movaps %s, %s" % (B4, C34))
    print("  dpps $0xf8, %s, %s" % (A3, C34))
    print("  blendps $0x07, %s, %s" % (C33, C34))
    C33.release()

    print("  addps %s, %s" % (C34, C3))
    C34.release()
    A3.release()

    A4 = SSERegister(log, "A4")
    issue_load(row_major_index(3, 0, 4)*4+offset(i, k, options.N)*16*4+spammOffsets.offset_block_dense, "A", A4)

    print("")
    print("  # Calculate C(4,:).")

    C41 = SSERegister(log, "C41")
    print("  movaps %s, %s" % (B1, C41))
    print("  dpps $0xf1, %s, %s" % (A4, C41))

    C42 = SSERegister(log, "C42")
    print("  movaps %s, %s" % (B2, C42))
    print("  dpps $0xf2, %s, %s" % (A4, C42))
    print("  blendps $0x01, %s, %s" % (C41, C42))
    C41.release()

    C43 = SSERegister(log, "C43")
    print("  movaps %s, %s" % (B3, C43))
    print("  dpps $0xf4, %s, %s" % (A4, C43))
    print("  blendps $0x03, %s, %s" % (C42, C43))
    C42.release()

    C44 = SSERegister(log, "C44")
    print("  movaps %s, %s" % (B4, C44))
    print("  dpps $0xf8, %s, %s" % (A4, C44))
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

  # We are done.
  return number_deactivated_products

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
  print("  addps 0x%x(C), %s" % (row_major_index(0, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense, C1))
  print("  addps 0x%x(C), %s" % (row_major_index(1, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense, C2))
  print("  addps 0x%x(C), %s" % (row_major_index(2, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense, C3))
  print("  addps 0x%x(C), %s" % (row_major_index(3, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense, C4))

  print("")
  print("  # Write out C(%d,%d) submatrix block." % (i+1, j+1))
  if options.no_store:
    print("  # skipped (command line option --no-store).")
    print("  #")
    print("  # movaps %s, 0x%x(C)" % (C1, row_major_index(0, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense))
    print("  # movaps %s, 0x%x(C)" % (C2, row_major_index(1, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense))
    print("  # movaps %s, 0x%x(C)" % (C3, row_major_index(2, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense))
    print("  # movaps %s, 0x%x(C)" % (C4, row_major_index(3, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense))
  else:
    print("  movaps %s, 0x%x(C)" % (C1, row_major_index(0, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense))
    print("  movaps %s, 0x%x(C)" % (C2, row_major_index(1, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense))
    print("  movaps %s, 0x%x(C)" % (C3, row_major_index(2, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense))
    print("  movaps %s, 0x%x(C)" % (C4, row_major_index(3, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense))

    if not options.hierarchical:
      last_store_offset.append(row_major_index(0, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense)
      last_store_offset.append(row_major_index(1, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense)
      last_store_offset.append(row_major_index(2, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense)
      last_store_offset.append(row_major_index(3, 0, 4)*4+offset(i, j, options.N)*16*4+spammOffsets.offset_block_dense)

  C1.release()
  C2.release()
  C3.release()
  C4.release()

# The main program.
parser = argparse.ArgumentParser(description =
"""This script generates a stream element kernel operating on 4x4 matrix
blocks. The kernel generated is written using assembly instructions assuming a
processor with the SSE instruction set.""")

parser.add_argument("-N",
    metavar = "N",
    help = "generate kernel for NxN matrix of 4x4 matrix blocks [default: %(default)s]",
    dest = "N",
    type = int,
    default = 1)

parser.add_argument("--name",
    metavar = "func",
    help = "set function name to \"func\" [default: %(default)s]",
    dest = "functionName",
    default = "stream_kernel")

parser.add_argument("--no-checks",
    help = "generate code without any norm checks [default: %(default)s]",
    action = "store_false",
    default = True,
    dest = "generate_checks")

parser.add_argument("--no-store",
    help = "do not generate store instructions (store C) [default: %(default)s]",
    action = "store_true",
    default = False)

parser.add_argument("--no-alpha",
    help = "generate a kernel for the special case of alpha = 1 [default: %(default)s]",
    action = "store_true",
    default = False,
    dest = "alphaOne")

parser.add_argument("--debug",
    help = "print out a lot of debugging information [default: %(default)s]",
    action = "store_true",
    default = False)

parser.add_argument("--deactivate-products",
    metavar = "D",
    help = "deactivate the first D products of the N^3 possible ones",
    type = int,
    default = 0)

parser.add_argument("--SSE",
    help = "set the SSE level from [1, 3, 4.1] [default: %(default)s]",
    metavar = "level",
    type = float,
    default = 1)

parser.add_argument("--Z-curve",
    help = """layout the multiply along a Z-curve as opposed to regular
row-major ordering [default: %(default)s]""",
    action = "store_true",
    default = False,
    dest = "Z_curve_ordering")

parser.add_argument("--hierarchical",
    help = "create a hierarchical kernel [default: %(default)s]",
    action = "store_true",
    default = False)

options = parser.parse_args()

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

if options.hierarchical:
  if options.N != 4:
    log.warning("Z-curve ordering only works with N = 4")
    options.N = 4

# Check SSE level.
if not options.SSE in [1, 3, 4.1]:
  log.error("unknown SSE level")
  sys.exit(1)

# Set the number of deactivated blocks.
if options.deactivate_products > options.N**3:
  log.error("there are only %d products in this kernel" % (options.N**3))
  sys.exit(1)

number_deactivated_products = options.deactivate_products

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
  alpha = SSERegister(log, "alpha", "%xmm0")
  print("#define %s %s" % (alpha.name, alpha.register))

# Get tolerance.
tolerance = SSERegister(log, "tolerance", "%xmm1")

print("#define %s %s" % (tolerance.name, tolerance.register))
print("#define multiply_stream %rsi")

print("")
print("# Define loop variables.")
print("#define index        %rax")
print("#define base_pointer %rdx")

print("")
print("# Define pointers to matrix data nodes in stream.")
print("#define A %r8")
print("#define B %r9")
print("#define C %r10")

if options.hierarchical:
  print("")
  print("# Define jump table variables.")
  print("#define jump_index         %r11")
  print("#define jump_index_temp    %r12")
  print("#define jump_index_base    %rcx")
  print("#define jump_index_base_32 %ecx")
  print("#define old_stack          %r13")

# Start the function prolog.
print("")
print("  # Function prolog.")
print("  .text")
print("  .balign 256")
print("  .global %s" % (options.functionName))
print("  .type %s, @function" % (options.functionName))

if options.hierarchical:
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
if options.hierarchical:
  print("  push jump_index")
  print("  push jump_index_temp")
  print("  push jump_index_base")
  print("  push old_stack")
print("  push A")
print("  push B")
print("  push C")

if options.hierarchical:
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

  norm_mask = SSERegister(log, "norm_mask")
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
print("  # Set the base pointer using sizeof(multiply_stream_t) = %d (0x%x)." % (spammOffsets.sizeof_multiply_stream_t, spammOffsets.sizeof_multiply_stream_t))
print("  imul $0x%x, base_pointer, base_pointer" % (spammOffsets.sizeof_multiply_stream_t))
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

# Initialize some other global variables.
last_store_offset = []

if options.hierarchical:
  print("")
  print("  # First level of hierarchy. Do some norm products [ A11*B11, A12*B21, A11*B12, A12*B22 ].")
  norm_1 = SSERegister(log, "norm_1")
  print("  movaps 0x%x(A), %s" % (spammOffsets.offset_norm_upper, norm_1))
  print("  mulps 0x%x(B), %s" % (spammOffsets.offset_norm_upper_transpose, norm_1))
  print("  # (normA*normB <= tolerance ? -1 : 0)")
  print("  cmpps $0x02, %s, %s # norm product <= tolerance?" % (tolerance, norm_1))
  print("  pshufb %s, %s" % (norm_mask, norm_1))
  print("  pmovmskb %s, jump_index" % (norm_1))
  print("")
  print("  # First level of hierarchy. Do some norm products [ A21*B11, A22*B21, A21*B12, A22*B22 ].")
  norm_2 = SSERegister(log, "norm_2")
  print("  movaps 0x%x(A), %s" % (spammOffsets.offset_norm_upper+4*4, norm_2))
  print("  mulps 0x%x(B), %s" % (spammOffsets.offset_norm_upper_transpose+4*4, norm_2))
  print("  # (normA*normB <= tolerance ? -1 : 0)")
  print("  cmpps $0x02, %s, %s # norm product <= tolerance?" % (tolerance, norm_2))
  print("  pshufb %s, %s" % (norm_mask, norm_2))
  print("  pmovmskb %s, jump_index_temp" % (norm_2))
  print("")
  print("  shl $4, jump_index_temp")
  print("  or jump_index_temp, jump_index")
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
  print("  # jump_index[4] = A21*B11")
  print("  # jump_index[5] = A22*B21")
  print("  # jump_index[6] = A21*B12")
  print("  # jump_index[7] = A22*B22")
  print("  #")
  print("  # where the array index i indicates 2^i.")

  print("")
  print("  # Jump table for first tier.")
  print("#if defined(__PIC__) || defined(__pic__)")
  print("  # For PIC (Position Independent Code) we need to perform some magic")
  print("  # when it comes to figuring out the relative addresses of the jump")
  print("  # table. Fortunatley on x86_64 we can use %rdi.")
  print("")
  print("  lea jump_table(%rip), jump_index_base")
  print("  mov (jump_index_base, jump_index, 4), jump_index_base_32")
  print("  movslq jump_index_base_32, jump_index_base # Expand pointer to 64-bit.")
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
      print("")
      print("  # Perform A_{11}*B_{11} = C_{11} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{11}*B_{11} + A_{12}*B_{21} = C_{11}")
      print("  # A_{11}*B_{12} + A_{12}*B_{22} = C_{12}")
      print("  # A_{21}*B_{11} + A_{22}*B_{21} = C_{21}")
      print("  # A_{21}*B_{12} + A_{22}*B_{22} = C_{22}")
      print("")
      print("  # Perform A_{12}*B_{21} = C_{11} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{13}*B_{31} + A_{14}*B_{41} = C_{11}")
      print("  # A_{13}*B_{32} + A_{14}*B_{42} = C_{12}")
      print("  # A_{23}*B_{31} + A_{24}*B_{41} = C_{21}")
      print("  # A_{23}*B_{32} + A_{24}*B_{42} = C_{22}")

      clearC(1, 1)
      number_deactivated_products = block_product(1, 1, 1, number_deactivated_products)
      number_deactivated_products = block_product(1, 2, 1, number_deactivated_products)
      number_deactivated_products = block_product(1, 3, 1, number_deactivated_products)
      number_deactivated_products = block_product(1, 4, 1, number_deactivated_products)
      writeC(1, 1)

      clearC(1, 2)
      number_deactivated_products = block_product(1, 1, 2, number_deactivated_products)
      number_deactivated_products = block_product(1, 2, 2, number_deactivated_products)
      number_deactivated_products = block_product(1, 3, 2, number_deactivated_products)
      number_deactivated_products = block_product(1, 4, 2, number_deactivated_products)
      writeC(1, 2)

      clearC(2, 1)
      number_deactivated_products = block_product(2, 1, 1, number_deactivated_products)
      number_deactivated_products = block_product(2, 2, 1, number_deactivated_products)
      number_deactivated_products = block_product(2, 3, 1, number_deactivated_products)
      number_deactivated_products = block_product(2, 4, 1, number_deactivated_products)
      writeC(2, 1)

      clearC(2, 2)
      number_deactivated_products = block_product(2, 1, 2, number_deactivated_products)
      number_deactivated_products = block_product(2, 2, 2, number_deactivated_products)
      number_deactivated_products = block_product(2, 3, 2, number_deactivated_products)
      number_deactivated_products = block_product(2, 4, 2, number_deactivated_products)
      writeC(2, 2)

    else:
      if (jump_index & 0x1) == 0:
        print("")
        print("  # Perform A_{11}*B_{11} = C_{11} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{11}*B_{11} + A_{12}*B_{21} = C_{11}")
        print("  # A_{11}*B_{12} + A_{12}*B_{22} = C_{12}")
        print("  # A_{21}*B_{11} + A_{22}*B_{21} = C_{21}")
        print("  # A_{21}*B_{12} + A_{22}*B_{22} = C_{22}")

        clearC(1, 1)
        number_deactivated_products = block_product(1, 1, 1, number_deactivated_products)
        number_deactivated_products = block_product(1, 2, 1, number_deactivated_products)
        writeC(1, 1)

        clearC(1, 2)
        number_deactivated_products = block_product(1, 1, 2, number_deactivated_products)
        number_deactivated_products = block_product(1, 2, 2, number_deactivated_products)
        writeC(1, 2)

        clearC(2, 1)
        number_deactivated_products = block_product(2, 1, 1, number_deactivated_products)
        number_deactivated_products = block_product(2, 2, 1, number_deactivated_products)
        writeC(2, 1)

        clearC(2, 2)
        number_deactivated_products = block_product(2, 1, 2, number_deactivated_products)
        number_deactivated_products = block_product(2, 2, 2, number_deactivated_products)
        writeC(2, 2)

      if (jump_index & 0x2) == 0:
        print("")
        print("  # Perform A_{12}*B_{21} = C_{11} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{13}*B_{31} + A_{14}*B_{41} = C_{11}")
        print("  # A_{13}*B_{32} + A_{14}*B_{42} = C_{12}")
        print("  # A_{23}*B_{31} + A_{24}*B_{41} = C_{21}")
        print("  # A_{23}*B_{32} + A_{24}*B_{42} = C_{22}")

        clearC(1, 1)
        number_deactivated_products = block_product(1, 3, 1, number_deactivated_products)
        number_deactivated_products = block_product(1, 4, 1, number_deactivated_products)
        writeC(1, 1)

        clearC(1, 2)
        number_deactivated_products = block_product(1, 3, 2, number_deactivated_products)
        number_deactivated_products = block_product(1, 4, 2, number_deactivated_products)
        writeC(1, 2)

        clearC(2, 1)
        number_deactivated_products = block_product(2, 3, 1, number_deactivated_products)
        number_deactivated_products = block_product(2, 4, 1, number_deactivated_products)
        writeC(2, 1)

        clearC(2, 2)
        number_deactivated_products = block_product(2, 3, 2, number_deactivated_products)
        number_deactivated_products = block_product(2, 4, 2, number_deactivated_products)
        writeC(2, 2)

    if (jump_index & (0x4 | 0x8)) == 0:
      print("")
      print("  # Perform A_{11}*B_{12} = C_{12} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{11}*B_{13} + A_{12}*B_{23} = C_{13}")
      print("  # A_{11}*B_{14} + A_{12}*B_{24} = C_{14}")
      print("  # A_{21}*B_{13} + A_{22}*B_{23} = C_{23}")
      print("  # A_{21}*B_{14} + A_{22}*B_{24} = C_{24}")
      print("")
      print("  # Perform A_{12}*B_{22} = C_{12} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{13}*B_{33} + A_{14}*B_{43} = C_{13}")
      print("  # A_{13}*B_{34} + A_{14}*B_{44} = C_{14}")
      print("  # A_{23}*B_{33} + A_{24}*B_{43} = C_{23}")
      print("  # A_{23}*B_{34} + A_{24}*B_{44} = C_{24}")

      clearC(1, 3)
      number_deactivated_products = block_product(1, 1, 3, number_deactivated_products)
      number_deactivated_products = block_product(1, 2, 3, number_deactivated_products)
      number_deactivated_products = block_product(1, 3, 3, number_deactivated_products)
      number_deactivated_products = block_product(1, 4, 3, number_deactivated_products)
      writeC(1, 3)

      clearC(1, 4)
      number_deactivated_products = block_product(1, 1, 4, number_deactivated_products)
      number_deactivated_products = block_product(1, 2, 4, number_deactivated_products)
      number_deactivated_products = block_product(1, 3, 4, number_deactivated_products)
      number_deactivated_products = block_product(1, 4, 4, number_deactivated_products)
      writeC(1, 4)

      clearC(2, 3)
      number_deactivated_products = block_product(2, 1, 3, number_deactivated_products)
      number_deactivated_products = block_product(2, 2, 3, number_deactivated_products)
      number_deactivated_products = block_product(2, 3, 3, number_deactivated_products)
      number_deactivated_products = block_product(2, 4, 3, number_deactivated_products)
      writeC(2, 3)

      clearC(2, 4)
      number_deactivated_products = block_product(2, 1, 4, number_deactivated_products)
      number_deactivated_products = block_product(2, 2, 4, number_deactivated_products)
      number_deactivated_products = block_product(2, 3, 4, number_deactivated_products)
      number_deactivated_products = block_product(2, 4, 4, number_deactivated_products)
      writeC(2, 4)

    else:
      if (jump_index & 0x4) == 0:
        print("")
        print("  # Perform A_{11}*B_{12} = C_{12} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{11}*B_{13} + A_{12}*B_{23} = C_{13}")
        print("  # A_{11}*B_{14} + A_{12}*B_{24} = C_{14}")
        print("  # A_{21}*B_{13} + A_{22}*B_{23} = C_{23}")
        print("  # A_{21}*B_{14} + A_{22}*B_{24} = C_{24}")

        clearC(1, 3)
        number_deactivated_products = block_product(1, 1, 3, number_deactivated_products)
        number_deactivated_products = block_product(1, 2, 3, number_deactivated_products)
        writeC(1, 3)

        clearC(1, 4)
        number_deactivated_products = block_product(1, 1, 4, number_deactivated_products)
        number_deactivated_products = block_product(1, 2, 4, number_deactivated_products)
        writeC(1, 4)

        clearC(2, 3)
        number_deactivated_products = block_product(2, 1, 3, number_deactivated_products)
        number_deactivated_products = block_product(2, 2, 3, number_deactivated_products)
        writeC(2, 3)

        clearC(2, 4)
        number_deactivated_products = block_product(2, 1, 4, number_deactivated_products)
        number_deactivated_products = block_product(2, 2, 4, number_deactivated_products)
        writeC(2, 4)

      if (jump_index & 0x8) == 0:
        print("")
        print("  # Perform A_{12}*B_{22} = C_{12} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{13}*B_{33} + A_{14}*B_{43} = C_{13}")
        print("  # A_{13}*B_{34} + A_{14}*B_{44} = C_{14}")
        print("  # A_{23}*B_{33} + A_{24}*B_{43} = C_{23}")
        print("  # A_{23}*B_{34} + A_{24}*B_{44} = C_{24}")

        clearC(1, 3)
        number_deactivated_products = block_product(1, 3, 3, number_deactivated_products)
        number_deactivated_products = block_product(1, 4, 3, number_deactivated_products)
        writeC(1, 3)

        clearC(1, 4)
        number_deactivated_products = block_product(1, 3, 4, number_deactivated_products)
        number_deactivated_products = block_product(1, 4, 4, number_deactivated_products)
        writeC(1, 4)

        clearC(2, 3)
        number_deactivated_products = block_product(2, 3, 3, number_deactivated_products)
        number_deactivated_products = block_product(2, 4, 3, number_deactivated_products)
        writeC(2, 3)

        clearC(2, 4)
        number_deactivated_products = block_product(2, 3, 4, number_deactivated_products)
        number_deactivated_products = block_product(2, 4, 4, number_deactivated_products)
        writeC(2, 4)

    if (jump_index & (0x10 | 0x20)) == 0:
      print("")
      print("  # Perform A_{21}*B_{11} = C_{21} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{31}*B_{11} + A_{32}*B_{21} = C_{31}")
      print("  # A_{31}*B_{12} + A_{32}*B_{22} = C_{32}")
      print("  # A_{41}*B_{11} + A_{42}*B_{21} = C_{41}")
      print("  # A_{41}*B_{12} + A_{42}*B_{22} = C_{42}")
      print("")
      print("  # Perform A_{22}*B_{21} = C_{21} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{33}*B_{31} + A_{34}*B_{41} = C_{31}")
      print("  # A_{33}*B_{32} + A_{34}*B_{42} = C_{32}")
      print("  # A_{43}*B_{31} + A_{44}*B_{41} = C_{41}")
      print("  # A_{43}*B_{32} + A_{44}*B_{42} = C_{42}")

      clearC(3, 1)
      number_deactivated_products = block_product(3, 1, 1, number_deactivated_products)
      number_deactivated_products = block_product(3, 2, 1, number_deactivated_products)
      number_deactivated_products = block_product(3, 3, 1, number_deactivated_products)
      number_deactivated_products = block_product(3, 4, 1, number_deactivated_products)
      writeC(3, 1)

      clearC(3, 2)
      number_deactivated_products = block_product(3, 1, 2, number_deactivated_products)
      number_deactivated_products = block_product(3, 2, 2, number_deactivated_products)
      number_deactivated_products = block_product(3, 3, 2, number_deactivated_products)
      number_deactivated_products = block_product(3, 4, 2, number_deactivated_products)
      writeC(3, 2)

      clearC(4, 1)
      number_deactivated_products = block_product(4, 1, 1, number_deactivated_products)
      number_deactivated_products = block_product(4, 2, 1, number_deactivated_products)
      number_deactivated_products = block_product(4, 3, 1, number_deactivated_products)
      number_deactivated_products = block_product(4, 4, 1, number_deactivated_products)
      writeC(4, 1)

      clearC(4, 2)
      number_deactivated_products = block_product(4, 1, 2, number_deactivated_products)
      number_deactivated_products = block_product(4, 2, 2, number_deactivated_products)
      number_deactivated_products = block_product(4, 3, 2, number_deactivated_products)
      number_deactivated_products = block_product(4, 4, 2, number_deactivated_products)
      writeC(4, 2)

    else:
      if (jump_index & 0x10) == 0:
        print("")
        print("  # Perform A_{21}*B_{11} = C_{21} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{31}*B_{11} + A_{32}*B_{21} = C_{31}")
        print("  # A_{31}*B_{12} + A_{32}*B_{22} = C_{32}")
        print("  # A_{41}*B_{11} + A_{42}*B_{21} = C_{41}")
        print("  # A_{41}*B_{12} + A_{42}*B_{22} = C_{42}")

        clearC(3, 1)
        number_deactivated_products = block_product(3, 1, 1, number_deactivated_products)
        number_deactivated_products = block_product(3, 2, 1, number_deactivated_products)
        writeC(3, 1)

        clearC(3, 2)
        number_deactivated_products = block_product(3, 1, 2, number_deactivated_products)
        number_deactivated_products = block_product(3, 2, 2, number_deactivated_products)
        writeC(3, 2)

        clearC(4, 1)
        number_deactivated_products = block_product(4, 1, 1, number_deactivated_products)
        number_deactivated_products = block_product(4, 2, 1, number_deactivated_products)
        writeC(4, 1)

        clearC(4, 2)
        number_deactivated_products = block_product(4, 1, 2, number_deactivated_products)
        number_deactivated_products = block_product(4, 2, 2, number_deactivated_products)
        writeC(4, 2)

      if (jump_index & 0x20) == 0:
        print("")
        print("  # Perform A_{22}*B_{21} = C_{21} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{33}*B_{31} + A_{34}*B_{41} = C_{31}")
        print("  # A_{33}*B_{32} + A_{34}*B_{42} = C_{32}")
        print("  # A_{43}*B_{31} + A_{44}*B_{41} = C_{41}")
        print("  # A_{43}*B_{32} + A_{44}*B_{42} = C_{42}")

        clearC(3, 1)
        number_deactivated_products = block_product(3, 3, 1, number_deactivated_products)
        number_deactivated_products = block_product(3, 4, 1, number_deactivated_products)
        writeC(3, 1)

        clearC(3, 2)
        number_deactivated_products = block_product(3, 3, 2, number_deactivated_products)
        number_deactivated_products = block_product(3, 4, 2, number_deactivated_products)
        writeC(3, 2)

        clearC(4, 1)
        number_deactivated_products = block_product(4, 3, 1, number_deactivated_products)
        number_deactivated_products = block_product(4, 4, 1, number_deactivated_products)
        writeC(4, 1)

        clearC(4, 2)
        number_deactivated_products = block_product(4, 3, 2, number_deactivated_products)
        number_deactivated_products = block_product(4, 4, 2, number_deactivated_products)
        writeC(4, 2)

    if (jump_index & (0x40 | 0x80)) == 0:
      print("")
      print("  # Perform A_{21}*B_{12} = C_{22} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{31}*B_{13} + A_{32}*B_{23} = C_{33}")
      print("  # A_{31}*B_{14} + A_{32}*B_{24} = C_{34}")
      print("  # A_{41}*B_{13} + A_{42}*B_{23} = C_{43}")
      print("  # A_{41}*B_{14} + A_{42}*B_{24} = C_{44}")
      print("")
      print("  # Perform A_{22}*B_{22} = C_{22} product with norm checks.")
      print("  #")
      print("  # On the basic matrix block level, this translates into:")
      print("  # A_{33}*B_{33} + A_{34}*B_{43} = C_{33}")
      print("  # A_{33}*B_{34} + A_{34}*B_{44} = C_{34}")
      print("  # A_{43}*B_{33} + A_{44}*B_{43} = C_{43}")
      print("  # A_{43}*B_{34} + A_{44}*B_{44} = C_{44}")

      clearC(3, 3)
      number_deactivated_products = block_product(3, 1, 3, number_deactivated_products)
      number_deactivated_products = block_product(3, 2, 3, number_deactivated_products)
      number_deactivated_products = block_product(3, 3, 3, number_deactivated_products)
      number_deactivated_products = block_product(3, 4, 3, number_deactivated_products)
      writeC(3, 3)

      clearC(3, 4)
      number_deactivated_products = block_product(3, 1, 4, number_deactivated_products)
      number_deactivated_products = block_product(3, 2, 4, number_deactivated_products)
      number_deactivated_products = block_product(3, 3, 4, number_deactivated_products)
      number_deactivated_products = block_product(3, 4, 4, number_deactivated_products)
      writeC(3, 4)

      clearC(4, 3)
      number_deactivated_products = block_product(4, 1, 3, number_deactivated_products)
      number_deactivated_products = block_product(4, 2, 3, number_deactivated_products)
      number_deactivated_products = block_product(4, 3, 3, number_deactivated_products)
      number_deactivated_products = block_product(4, 4, 3, number_deactivated_products)
      writeC(4, 3)

      clearC(4, 4)
      number_deactivated_products = block_product(4, 1, 4, number_deactivated_products)
      number_deactivated_products = block_product(4, 2, 4, number_deactivated_products)
      number_deactivated_products = block_product(4, 3, 4, number_deactivated_products)
      number_deactivated_products = block_product(4, 4, 4, number_deactivated_products)
      writeC(4, 4)

    else:
      if (jump_index & 0x40) == 0:
        print("")
        print("  # Perform A_{21}*B_{12} = C_{22} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{31}*B_{13} + A_{32}*B_{23} = C_{33}")
        print("  # A_{31}*B_{14} + A_{32}*B_{24} = C_{34}")
        print("  # A_{41}*B_{13} + A_{42}*B_{23} = C_{43}")
        print("  # A_{41}*B_{14} + A_{42}*B_{24} = C_{44}")

        clearC(3, 3)
        number_deactivated_products = block_product(3, 1, 3, number_deactivated_products)
        number_deactivated_products = block_product(3, 2, 3, number_deactivated_products)
        writeC(3, 3)

        clearC(3, 4)
        number_deactivated_products = block_product(3, 1, 4, number_deactivated_products)
        number_deactivated_products = block_product(3, 2, 4, number_deactivated_products)
        writeC(3, 4)

        clearC(4, 3)
        number_deactivated_products = block_product(4, 1, 3, number_deactivated_products)
        number_deactivated_products = block_product(4, 2, 3, number_deactivated_products)
        writeC(4, 3)

        clearC(4, 4)
        number_deactivated_products = block_product(4, 1, 4, number_deactivated_products)
        number_deactivated_products = block_product(4, 2, 4, number_deactivated_products)
        writeC(4, 4)

      if (jump_index & 0x80) == 0:
        print("")
        print("  # Perform A_{22}*B_{22} = C_{22} product with norm checks.")
        print("  #")
        print("  # On the basic matrix block level, this translates into:")
        print("  # A_{33}*B_{33} + A_{34}*B_{43} = C_{33}")
        print("  # A_{33}*B_{34} + A_{34}*B_{44} = C_{34}")
        print("  # A_{43}*B_{33} + A_{44}*B_{43} = C_{43}")
        print("  # A_{43}*B_{34} + A_{44}*B_{44} = C_{44}")

        clearC(3, 3)
        number_deactivated_products = block_product(3, 3, 3, number_deactivated_products)
        number_deactivated_products = block_product(3, 4, 3, number_deactivated_products)
        writeC(3, 3)

        clearC(3, 4)
        number_deactivated_products = block_product(3, 3, 4, number_deactivated_products)
        number_deactivated_products = block_product(3, 4, 4, number_deactivated_products)
        writeC(3, 4)

        clearC(4, 3)
        number_deactivated_products = block_product(4, 3, 3, number_deactivated_products)
        number_deactivated_products = block_product(4, 4, 3, number_deactivated_products)
        writeC(4, 3)

        clearC(4, 4)
        number_deactivated_products = block_product(4, 3, 4, number_deactivated_products)
        number_deactivated_products = block_product(4, 4, 4, number_deactivated_products)
        writeC(4, 4)

    print("")
    print("  # Done.")
    print("  jmp loop_end")

else:
  for i in range(options.N):
    for j in range(options.N):
      clearC(i+1, j+1)
      for k in range(options.N):
        number_deactivated_products = block_product(i+1, k+1, j+1, number_deactivated_products)
      writeC(i+1, j+1)

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
if options.hierarchical:
  print("")
  print("  # Restore old stack.")
  print("  mov old_stack, %rsp")
print("")
print("  # Pop registers from stack.")
print("  pop C")
print("  pop B")
print("  pop A")
if options.hierarchical:
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

if options.hierarchical:
  norm_mask.release()

# Start function epilog.
print("")
print("  # Function epilog.")
print("  .size %s, .-%s" % (options.functionName, options.functionName))

# For debugging.
log.debug("remaining variables assigned as SSERegister: %s" % (SSERegister.variables))

if len(SSERegister.variables) > 0:
  log.error("still assigned variables as SSERegister: %s" % (SSERegister.variables))
  sys.exit(1)
