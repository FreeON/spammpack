AC_INIT([libchunk],
        [1],
        [nicolas.bock@freeon.org],
        [],
        [http://www.freeon.org/])
AM_INIT_AUTOMAKE([dist-bzip2 no-dist-gzip silent-rules])
AM_SILENT_RULES([yes])

# Find a C compiler.
AC_PROG_CC_STDC

# Find assembler.
AM_PROG_AS

# Turn on asserts.
AC_HEADER_ASSERT

# Initialize libtool.
LT_INIT

# Find OpenMP for C.
AC_OPENMP

# Check for getopt.
AC_CHECK_HEADER([getopt.h], [], [AC_MSG_FAILURE([need getopt.h])])
AC_SEARCH_LIBS([getopt_long], [],
  [], [AC_MSG_FAILURE([need getopt_long()])])

# Test for math.
AC_SEARCH_LIBS([ceil], [m], [], [AC_MSG_FAILURE([need the math library])])

# Test for clock_gettime.
AC_SEARCH_LIBS([clock_gettime], [rt],
  [],
  [AC_MSG_FAILURE([need clock_gettime()])])

# Check for supported clocks. Thanks to Jonathan Lifflander
# <jliffl2@illinois.edu> for suggesting this test.
for CLOCKTYPE in \
  CLOCK_MONOTONIC_RAW \
  CLOCK_MONOTONIC \
  CLOCK_REALTIME; do
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <time.h>]],
    [[struct timespec p; clock_gettime(${CLOCKTYPE}, &p);]])],
    [break], [CLOCKTYPE="unset"])
done
if test "${CLOCKTYPE}" = "unset"; then
  AC_MSG_FAILURE([no known clock type])
fi
AC_DEFINE_UNQUOTED([CLOCKTYPE], [${CLOCKTYPE}], [The clock type to use with clock_gettime()])

AC_ARG_ENABLE([mic],
  [AS_HELP_STRING([--enable-mic],
  [enable Xeon PHI])],
  [enable_mic=${enableval}],
  [enable_mic="no"])
if [[ ${enable_mic} = "yes" ]]; then
  AC_DEFINE([MIC_ALLOC], [1], [Enable huge page allocation for MIC])
  AC_SUBST([MIC_CFLAGS], [-mmic])
  AC_SUBST([MIC_LDFLAGS], [-static])
fi

AC_ARG_ENABLE([set-num-threads],
  [AS_HELP_STRING([--enable-set-num-threads],
  [set the number of threads in OpenMP parallel sections to N])],
  [enable_set_threads=${enableval}],
  [enable_set_threads="no"])
if [[ ${enable_set_threads} != "no" ]]; then
  AC_DEFINE_UNQUOTED([FORCE_OMP_NUM_THREADS], [${enable_set_threads}],
    [The number of threads in OpenMP parallel sections])
fi

AC_MSG_CHECKING([which chunk implementation to use])
AC_ARG_ENABLE([chunk-implementation],
  [AS_HELP_STRING([--enable-chunk-implementation],
  [choose the chunk implementation (tiled, tree), default = tree])],
  [enable_chunk_implementation=$enableval],
  [enable_chunk_implementation="tree"])
case ${enable_chunk_implementation} in
  "tiled" )
    AC_MSG_RESULT([tiled chunks])
    ;;
  "tree" )
    AC_MSG_RESULT([tree chunks])
    ;;
  * )
    AC_MSG_FAILURE([unknown chunk implementation: ${enable_chunk_implementation}])
    ;;
esac
AC_DEFINE_UNQUOTED([CHUNK_IMPLEMENTATION],
  [${enable_chunk_implementation}], [The chunk implementation])

AC_MSG_CHECKING([which block multiply implementation to use])
AC_ARG_ENABLE([block-multiply],
  [AS_HELP_STRING([--enable-block-multiply],
  [choose the dense multiply implementation (block, blas), default = block])],
  [enable_block_multiply=$enableval],
  [enable_block_multiply="block"])
case ${enable_block_multiply} in
  "block" )
    AC_MSG_RESULT([block])
    AC_DEFINE([BLOCK_MULTIPLY], [1], [Use simple nested loop matrix multiply.])
    ;;
  "blas" )
    AC_MSG_RESULT([blas])
    AC_DEFINE([BLOCK_BLAS], [1], [Use gemm() matrix multiply.])
    ;;
  *) AC_MSG_FAILURE([unknown multiply implementation]) ;;
esac

AC_ARG_ENABLE([lapack],
  [AS_HELP_STRING([--enable-lapack],
  [enable external linkage to blas/lapack])],
  [use_external_lapack=$enableval],
  [use_external_lapack="no"])

if [[ "${use_external_lapack}" = "yes" ]] || [[ "${enable_block_multiply}" = "blas" ]]; then
  # Find dgemm() for matrix block products.
  for DGEMM in dgemm dgemm_ dgemm__; do
    AC_SEARCH_LIBS([${DGEMM}], [blas], [break], [DGEMM="unset"])
  done
  if test "${DGEMM}" != "unset"; then
    AC_DEFINE_UNQUOTED([DGEMM], [${DGEMM}], [The call to BLAS dgemm() including
                        all necessary underscores.])
  fi

  # Find dsyev() for spectral bounds.
  for DSYEV in dsyev dsyev_ dsyev__; do
    AC_SEARCH_LIBS([${DSYEV}], [lapack], [break], [DSYEV="unset"])
  done
  if test "${DSYEV}" != "unset"; then
    AC_DEFINE_UNQUOTED([DSYEV], [${DSYEV}], [The call to dsyev() including all
                        necessary underscores.])
  fi
fi

AC_ARG_ENABLE([chunk-tree-max-tier],
  [AS_HELP_STRING([--enable-chunk-tree-max-tier],
  [set the maximum tier at which new tasks are created when recursing the
   chunk tree (default = 6)])],
  [enable_chunk_tree_max_tier=$enableval],
  [enable_chunk_tree_max_tier=6])
if [[ ${enable_chunk_tree_max_tier} -gt 0 ]]; then
  AC_DEFINE_UNQUOTED([CHUNK_TREE_MAX_TIER], [${enable_chunk_tree_max_tier}],
    [The maximum tier at which new tasks are created when recursing the chunk
     tree.])
fi

AC_ARG_ENABLE([chunk-block-transpose],
  [AS_HELP_STRING([--enable-chunk-block-transpose],
  [set the minimum size N of the block matrix that will transpose a matrix to
   allow for loop vectorization])],
  [enable_chunk_block_transpose=$enableval],
  [enable_chunk_block_transpose="32"])
if [[ ${enable_chunk_block_transpose} -gt 0 ]]; then
  AC_DEFINE_UNQUOTED([CHUNK_BLOCK_TRANSPOSE], [${enable_chunk_block_transpose}],
    [The minimum size N of the block matrix that will transpose a matrix to
     allow for loop vectorization.])
fi

AC_ARG_ENABLE([sleep],
  [AS_HELP_STRING([--enable-sleep],
  [add N nanoseconds of sleep during the block multiply])],
  [enable_sleep=$enableval],
  [enable_sleep=0])
if [[ ${enable_sleep} -gt 0 ]]; then
  AC_DEFINE_UNQUOTED([BLOCK_SLEEP], [${enable_sleep}], [Number of nanoseconds
    to sleep after block multiply.])
fi

AC_ARG_ENABLE([complexity],
  [AS_HELP_STRING([--enable-complexity],
  [count the product complexity in the tree multiply])],
  [enable_complexity=$enableval],
  [enable_complexity="no"])
if [[ ${enable_complexity} = "yes" ]]; then
  AC_DEFINE([MEASURE_COMPLEXITY], [1], [Measure the product complexity])
fi

# Enable debugging.
AC_MSG_CHECKING([whether to turn on full debugging output])
AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug],
  [enable a lot of extra output])],
  [enable_debug=$enableval],
  [enable_debug="no"])
if test "${enable_debug}" = "yes"; then
  AC_DEFINE([DEBUG_OUTPUT], [1], [Enable a lot of extra output])
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

# Enable matrix printing.
AC_MSG_CHECKING([whether to print matrices])
AC_ARG_ENABLE([debug-matrix],
  [AS_HELP_STRING([--enable-debug-matrix],
  [print the matrices for debugging])],
  [enable_debug=$enableval],
  [enable_debug="no"])
if test "${enable_debug}" = "yes"; then
  AC_DEFINE([DEBUG_OUTPUT_MATRIX], [1], [Print the matrices for debugging])
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
Makefile
src/Makefile
tests/Makefile
])

AC_OUTPUT
