AC_INIT([matmul],
        [6],
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

AC_SEARCH_LIBS([clock_gettime], [rt],
               [],
               [AC_MSG_FAILURE([need clock_gettime()])])

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
                  [Set the default load balancer.])],
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
                  [Set the absolute tolerance of the verification])],
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

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile Doxyfile])

AC_OUTPUT
