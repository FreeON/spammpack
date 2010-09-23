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

# Generate offsets.
print
for i in range(options.N):
  for j in range(options.N):
    print "#define A_OFFSET_%d%d (%d*%d+%d)*64 // %d" % (i+1, j+1, i, options.N, j, (i*options.N+j)*64)

print
for i in range(options.N):
  for j in range(options.N):
    print "#define B_OFFSET_%d%d (%d*%d+%d)*16 // %d" % (i+1, j+1, i, options.N, j, (i*options.N+j)*16)

print
for i in range(options.N):
  for j in range(options.N):
    print "#define C_OFFSET_%d%d (%d*%d+%d)*16 // %d" % (i+1, j+1, i, options.N, j, (i*options.N+j)*16)

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
print "  short int i;"
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

for i in range(options.N):
  for j in range(options.N):
    print
    sys.stdout.write("    if (")
    for k in range(options.N):
      norm_index_A = i*options.N+k
      norm_index_B = k*options.N+j+options.N**2
      sys.stdout.write("norm[%d]*norm[%d] >= tolerance" % (norm_index_A, norm_index_B))
      if k < options.N-1:
        sys.stdout.write(" &&\n")
        sys.stdout.write("        ")
      else:
        sys.stdout.write(")\n")

    print "    {"
    for k in range(options.N):
      print "      /* A(%d,%d)*B(%d,%d) = C(%d,%d). */" % (i+1, k+1, k+1, j+1, i+1, j+1)
      print "      for (i = 0; i < 4; i++)"
      print "      {"
      print "        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_%d%d]);" % (i+1, k+1)
      print "        B_row = _mm_load_ps(&B[0*4+B_OFFSET_%d%d]);" % (k+1, j+1)
      if k == 0:
        print "        C_row[i] = _mm_mul_ps(A_element, B_row);"
      else:
        print "        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);"
      print
      print "        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_%d%d]);" % (i+1, k+1)
      print "        B_row = _mm_load_ps(&B[1*4+B_OFFSET_%d%d]);" % (k+1, j+1)
      print "        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);"
      print
      print "        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_%d%d]);" % (i+1, k+1)
      print "        B_row = _mm_load_ps(&B[2*4+B_OFFSET_%d%d]);" % (k+1, j+1)
      print "        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);"
      print
      print "        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_%d%d]);" % (i+1, k+1)
      print "        B_row = _mm_load_ps(&B[3*4+B_OFFSET_%d%d]);" % (k+1, j+1)
      print "        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);"
      print "      }"
      print

    print "      /* Store C block. */"
    print "      for (i = 0; i < 4; i++)"
    print "      {"
    print "        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);"
    print "        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_%d%d]), C_row[i]);" % (i+1, j+1)
    print "        _mm_store_ps(&C[i*4+C_OFFSET_%d%d], C_row[i]);" % (i+1, j+1)
    print "      }"
    print "    }"

print "  }"
print "}"
