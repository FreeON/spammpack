#!/usr/bin/python
#
# Generate C SSE intrinsics code for a kernel operating on a 4x4 blocks.

import math, optparse, sys

parser = optparse.OptionParser(description =
"""This script generates a stream element kernel operating on 4x4 matrix
blocks. The kernel is written using C/C++ SSE intrinsics""")

parser.add_option("-N",
    metavar = "N",
    help = "generate kernel for NxN matrix of 4x4 matrix blocks",
    dest = "N",
    type = "int",
    default = 1)

parser.add_option("--name",
    metavar = "func",
    help = "set function name to \"func\"",
    dest = "functionName",
    type = "string",
    default = "stream_kernel")

( options, arguments ) = parser.parse_args()

# Check N.
d = int(math.log(options.N)/math.log(2))
if 2**d != options.N:
  print "N needs to be a power of 2"
  sys.exit(1)

# Generate C code.
print "#include <stdlib.h>"
print "#include <xmmintrin.h>"

print
print "struct multiply_stream_t"
print "{"
print "  float *A_block;"
print "  float *B_block;"
print "  float *C_block;"
print "  float  norm[%d];" % (2*options.N**2)
print "};"
print
print "void"
print "%s (const unsigned int number_stream_elements," % (options.functionName)
print "    float alpha,"
print "    float tolerance,"
print "    struct multiply_stream_t *multiply_stream)"
print "{"

# Generate offsets.
print
sys.stdout.write("  static const int A_OFFSET[%d][%d] =\n" % (options.N, options.N))
sys.stdout.write("  {\n")
for i in range(options.N):
  sys.stdout.write("    {\n")
  for j in range(options.N):
    sys.stdout.write("      %d" % ((i*options.N+j)*64))
    if j < options.N-1:
      sys.stdout.write(",")
    sys.stdout.write(" /* (%d*%d+%d)*64 */\n" % (i, options.N, j))
  sys.stdout.write("    }")
  if i < options.N-1:
    sys.stdout.write(",")
  sys.stdout.write("\n")
sys.stdout.write("  };\n")

print
sys.stdout.write("  static const int B_OFFSET[%d][%d] =\n" % (options.N, options.N))
sys.stdout.write("  {\n")
for i in range(options.N):
  sys.stdout.write("    {\n")
  for j in range(options.N):
    sys.stdout.write("      %d" % ((i*options.N+j)*16))
    if j < options.N-1:
      sys.stdout.write(",")
    sys.stdout.write(" /* (%d*%d+%d)*16 */\n" % (i, options.N, j))
  sys.stdout.write("    }")
  if i < options.N-1:
    sys.stdout.write(",")
  sys.stdout.write("\n")
sys.stdout.write("  };\n")

print
sys.stdout.write("  static const int C_OFFSET[%d][%d] =\n" % (options.N, options.N))
sys.stdout.write("  {\n")
for i in range(options.N):
  sys.stdout.write("    {\n")
  for j in range(options.N):
    sys.stdout.write("      %d" % ((i*options.N+j)*16))
    if j < options.N-1:
      sys.stdout.write(",")
    sys.stdout.write(" /* (%d*%d+%d)*16 */\n" % (i, options.N, j))
  sys.stdout.write("    }")
  if i < options.N-1:
    sys.stdout.write(",")
  sys.stdout.write("\n")
sys.stdout.write("  };\n")

print
print "  short int i, j, k, l;"
print "  unsigned int stream_index;"
print "  unsigned int max_stream_index;"
print
print "  float *restrict A;"
print "  float *restrict B;"
print "  float *restrict C;"
print
print "  float *restrict norm;"
print
print "  __m128 alpha_row;"
print
print "  __m128 A_element;"
print "  __m128 B_row;"
print "  __m128 C_row[4];"
print
print "  /* Divide number of stream elements by %d to simulate stride of %d. */" % (options.N**3, options.N**3)
print "  max_stream_index = number_stream_elements/%d;" % (options.N**3)
print
print "  alpha_row = _mm_set1_ps(alpha);"
print
print "  for (stream_index = 0; stream_index < max_stream_index; stream_index++)"
print "  {"
print "    /* Load pointers to matrix data blocks. */"
print "    A = multiply_stream[stream_index].A_block;"
print "    B = multiply_stream[stream_index].B_block;"
print "    C = multiply_stream[stream_index].C_block;"
print "    norm = multiply_stream[stream_index].norm;"

print
print "    for (i = 0; i < %d; i++) {" % (options.N)
print "      for (j = 0; j < %d; j++)" % (options.N)
print "      {"
print "        C_row[0] = _mm_setzero_ps();"
print "        C_row[1] = _mm_setzero_ps();"
print "        C_row[2] = _mm_setzero_ps();"
print "        C_row[3] = _mm_setzero_ps();"
print
print "        for (k = 0; k < %d; k++)" % (options.N)
print "        {"
print "          if (norm[i*%d+k]*norm[k*%d+j+%d] < tolerance) { continue; }" % (options.N, options.N, options.N**2)
print "          else"
print "          {"
print "            for (l = 0; l < 4; l++)"
print "            {"
print "              A_element = _mm_load_ps(&A[(l*4+0)*4+A_OFFSET[i][k]]);"
print "              B_row = _mm_load_ps(&B[0*4+B_OFFSET[k][j]]);"
print "              C_row[l] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[l]);"
print
print "              A_element = _mm_load_ps(&A[(l*4+1)*4+A_OFFSET[i][k]]);"
print "              B_row = _mm_load_ps(&B[1*4+B_OFFSET[k][j]]);"
print "              C_row[l] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[l]);"
print
print "              A_element = _mm_load_ps(&A[(l*4+2)*4+A_OFFSET[i][k]]);"
print "              B_row = _mm_load_ps(&B[2*4+B_OFFSET[k][j]]);"
print "              C_row[l] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[l]);"
print
print "              A_element = _mm_load_ps(&A[(l*4+3)*4+A_OFFSET[i][k]]);"
print "              B_row = _mm_load_ps(&B[3*4+B_OFFSET[k][j]]);"
print "              C_row[l] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[l]);"
print "            }"
print "          }"
print "        }"
print
print "        /* Store C block. */"
print "        for (l = 0; l < 4; l++)"
print "        {"
print "          C_row[l] = _mm_mul_ps(alpha_row, C_row[l]);"
print "          C_row[l] = _mm_add_ps(_mm_load_ps(&C[l*4+C_OFFSET[i][j]]), C_row[l]);"
print "          _mm_store_ps(&C[l*4+C_OFFSET[i][j]], C_row[l]);"
print "        }"
print "      }"
print "    }"
print "  }"
print "}"
