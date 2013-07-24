AC_INIT([matmul], [4], [nicolasbock@gmail.com])
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
  AC_SUBST([LOADBALANCER_FLAGS], [-module CommonLBs])
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

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile Doxyfile])

AC_OUTPUT
