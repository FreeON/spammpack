AC_INIT([matmul], [1], [nicolasbock@gmail.com])
AM_INIT_AUTOMAKE([silent-rules])
AM_SILENT_RULES([yes])

AC_PROG_CC

AC_ARG_VAR([CHARMC], [The Charm++ compiler wrapper])
AC_CHECK_PROG([CHARMC], [charmc], [charmc])

AC_PROG_CXX([${CHARMC}])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])

AC_OUTPUT
