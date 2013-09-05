AC_INIT([matmul],
        [15],
        [nicolas.bock@freeon.org],
        [matmul],
        [http://www.freeon.org/])
AM_INIT_AUTOMAKE([dist-bzip2 no-dist-gzip silent-rules])
AM_SILENT_RULES([yes])

# Find doxygen and related programs.
AC_ARG_VAR([DOXYGEN], [The doxygen executable (for building the documentation)])
AC_CHECK_PROG([DOXYGEN], [doxygen], [doxygen])

AC_CHECK_PROG([DOT], [dot], [dot])
if test -n ${DOT}; then
  have_dot="yes"
else
  have_dot="no"
fi
AC_SUBST([PROJECT_NAME], [${PACKAGE_NAME}])
AC_SUBST([PROJECT_NUMBER], [${PACKAGE_VERSION}])
AC_SUBST([HAVE_DOT], [${have_dot}])

# The doxygen configuration file is generated.
AC_SUBST([CONFIG_STATUS_DEPENDENCIES], ['$(top_srcdir)/Doxyfile.in'])

# Find a python interpreter.
AM_PATH_PYTHON([2.7], [], [:])
AM_CONDITIONAL([HAVE_PYTHON], [test "${PYTHON}" != ":"])

# Find the Charm++ compiler wrapper.
AC_ARG_VAR([CHARMC], [The Charm++ compiler wrapper])
if test -z "${CHARMC}"; then
  CHARMC="charmc"
fi
CXX=${CHARMC}

# Find the real compiler.
AC_PROG_CC
AC_PROG_CXX

# Compiler tests should be done using the real C compiler and not the charmc
# wrapper.
AC_LANG([C])

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

# Find dgemm() for matrix block products.
for DGEMM in dgemm dgemm_ dgemm__; do
  AC_SEARCH_LIBS([${DGEMM}], [blas], [break], [DGEMM="unset"])
done
if test "${DGEMM}" != "unset"; then
  AC_DEFINE_UNQUOTED([DGEMM], [${DGEMM}], [The call to BLAS dgemm() including
                      all necessary underscores.])
fi

# Enable load balancers.
AC_ARG_ENABLE([load-balance],
  [AS_HELP_STRING([--enable-load-balance],
  [enable a common set of load balancers during compilation])],
  [enable_LB=$enableval],
  [enable_LB="no"])
if test "${enable_LB}" = "yes"; then
  AC_SUBST([LOADBALANCER_FLAGS], ["-module CommonLBs"])
fi

# Enable the default load balancer.
AC_ARG_ENABLE([default-LB],
  [AS_HELP_STRING([--enable-default-LB],
  [set the default load balancer (default = GreedyCommLB)])],
  [default_LB=$enableval], [])
if test -n "${default_LB}"; then
  if test "${default_LB}" = "yes"; then
    default_LB="GreedyCommLB"
  fi
  AC_SUBST([DEFAULT_LB_FLAGS], ["-balancer ${default_LB}"])
fi

# Enable charmdebug.
AC_ARG_ENABLE([charmdebug],
  [AS_HELP_STRING([--enable-charmdebug],
  [enable the use of charmdebug on the executable])],
  [enable_charmdebug=$enableval],
  [enable_charmdebug="no"])
if test "${enable_charmdebug}" = "yes"; then
  AC_SUBST([CHARMDEBUG_FLAGS], ["-memory charmdebug"])
fi

# Enable debugging.
AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug],
  [enable a lot of extra output])],
  [enable_debug=$enableval],
  [enable_debug="no"])
if test "${enable_debug}" = "yes"; then
  AC_DEFINE([DEBUG_OUTPUT], [1], [Enable a lot of extra output])
fi

# Print some matrices for debugging.
AC_ARG_ENABLE([print-matrices],
  [AS_HELP_STRING([--enable-print-matrices],
  [print matrices for debugging])],
  [enable_print_matrices=$enableval],
  [enable_print_matrices="no"])
if test "${enable_print_matrices}" = "yes"; then
  AC_DEFINE([PRINT_MATRICES], [1], [Print matrices for debugging])
fi

# Enable pruning of the convolution. This is done by calling ckDestroy() on
# convolution elements.
AC_ARG_ENABLE([pruning],
  [AS_HELP_STRING([--enable-pruning],
  [enable pruning of convolution space])],
  [enable_pruning=$enableval],
  [enable_pruning="no"])
if test "${enable_pruning}" = "yes"; then
  AC_DEFINE([PRUNE_CONVOLUTION], [1], [Enable pruning of convolution space])
fi

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile Doxyfile])

AC_OUTPUT
