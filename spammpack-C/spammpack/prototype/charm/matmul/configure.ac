AC_INIT([matmul],
        [7],
        [nicolas.bock@freeon.org],
        [matmul],
        [http://www.freeon.org/])
AM_INIT_AUTOMAKE([dist-bzip2 no-dist-gzip silent-rules])
AM_SILENT_RULES([yes])

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

AC_SUBST([CONFIG_STATUS_DEPENDENCIES], ['$(top_srcdir)/Doxyfile.in'])

AC_ARG_VAR([CHARMC], [The Charm++ compiler wrapper])
AC_CHECK_PROG([CHARMC], [charmc], [charmc])

AC_PROG_CC
AC_PROG_CXX([${CHARMC}])

AC_LANG([C])

AC_SEARCH_LIBS([clock_gettime], [rt],
  [],
  [AC_MSG_FAILURE([need clock_gettime()])])

for CLOCKTYPE in CLOCK_MONOTONIC_RAW CLOCK_MONOTONIC CLOCK_REALTIME; do
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <time.h>]],
    [[struct timespec p; clock_gettime(${CLOCKTYPE}, &p);]])],
    [break], [CLOCKTYPE="unset"])
done
if test ${CLOCKTYPE} = "unset"; then
  AC_MSG_FAILURE([no known clock type])
fi
AC_DEFINE_UNQUOTED([CLOCKTYPE], [${CLOCKTYPE}], [The clock type to use with clock_gettime()])

AC_ARG_ENABLE([load-balance],
  [AS_HELP_STRING([--enable-load-balance],
                  [Enable a common set of load balancers during compile])],
  [enable_LB=$enableval],
  [enable_LB="no"])
if test "${enable_LB}" = "yes"; then
  AC_SUBST([LOADBALANCER_FLAGS], ["-module CommonLBs"])
fi

AC_ARG_WITH([default-LB],
  [AS_HELP_STRING([--with-default-LB],
  [Set the default load balancer (default = GreedyCommLB).])],
  [default_LB=$withval], [])
if test -n "${default_LB}"; then
  if test "${default_LB}" = "yes"; then
    default_LB="GreedyCommLB"
  fi
  AC_SUBST([DEFAULT_LB_FLAGS], ["-balancer ${default_LB}"])
fi

AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug],
  [Enable a lot of extra output])],
  [enable_debug=$enableval],
  [enable_debug="no"])
if test "${enable_debug}" = "yes"; then
  AC_DEFINE([DEBUG_OUTPUT], [1], [Enable a lot of extra output])
fi

AC_ARG_ENABLE([verify],
  [AS_HELP_STRING([--enable-verify],
  [Verify multiplication result])],
  [enable_verify=$enableval],
  [enable_verify="no"])
if test "${enable_verify}" = "yes"; then
  AC_DEFINE([VERIFY_MULTIPLY], [1], [Verify multiplication result])
fi

AC_ARG_WITH([verify-tolerance],
  [AS_HELP_STRING([--with-verify-tolerance],
  [Set the absolute tolerance of the verification (default = 1.0e-10).])],
  [verify_tolerance=$withval],
  [verify_tolerance="1.0e-10"])
AC_DEFINE_UNQUOTED([VERIFY_TOLERANCE], [${verify_tolerance}], [The absoluve
                    toleranc of the multiply verification])

AC_ARG_ENABLE([reduction-test],
  [AS_HELP_STRING([--enable-reduction-test],
  [Instead of multiplying, run a simple reduction test.])],
  [reduction_test=$enableval],
  [reduction_test="no"])
if test "${reduction_test}" = "yes"; then
  AC_DEFINE([REDUCTION_TEST], [1], [Instead of multiplying, run a simple reduction test.])
fi

AC_ARG_ENABLE([direct-multiply],
  [AS_HELP_STRING([--enable-direct-multiply],
  [Construct the convolution directly in the main program.])],
  [direct_multiply=$enableval],
  [direct_multiply="no"])
if test "${direct_multiply}" = "yes"; then
  AC_DEFINE([DIRECT_MULTIPLY], [1], [Construct the convolution curve directly in the main program.])
fi

AC_ARG_ENABLE([reduction-target],
  [AS_HELP_STRING([--enable-reduction-target],
  [Use an explicit reduction target in the multiply instead of suspending the
   calling thread.])],
  [reduction_target=$enableval],
  [reduction_target="no"])
if test "${reduction_target}" = "yes"; then
  AC_DEFINE([USE_REDUCTION_TARGET], [1], [Use an explicit reduction target in
             the multiply instead of suspending the calling thread.])
fi

AC_ARG_WITH([decay],
  [AS_HELP_STRING([--with-decay],
  [Set the decay constant for matrics with decay (default = 0.1).])],
  [matrix_decay=$withval],
  [matrix_decay="0.1"])
AC_DEFINE_UNQUOTED([MATRIX_DECAY], [${matrix_decay}], [Set the decay constant
                    for matrics with decay.])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile Doxyfile])

AC_OUTPUT
