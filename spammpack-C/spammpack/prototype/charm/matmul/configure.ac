AC_INIT([matmul], [1], [nicolasbock@gmail.com])
AM_INIT_AUTOMAKE([silent-rules])
AM_SILENT_RULES([yes])

AC_PROG_CC

AC_ARG_VAR([CHARMC], [The Charm++ compiler wrapper])
AC_CHECK_PROG([CHARMC], [charmc], [charmc])
AC_PROG_CXX([${CHARMC}])

AC_ARG_ENABLE([load-balance],
  [AS_HELP_STRING([--enable-load-balance],
                  [Enable a common set of load balancers during compile])],
  [enable_LB=$enableval],
  [enable_LB="no"])
if test "${enable_LB}" = "yes"; then
  AC_SUBST([LOADBALANCER_FLAGS], [-module CommonLBs])
fi

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])

AC_OUTPUT
