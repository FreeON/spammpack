AC_INIT([SpAMM], [0.1], [nbock@lanl.gov])
AM_INIT_AUTOMAKE

AC_CONFIG_MACRO_DIR([m4])

dnl Identify compilers.
AC_PROG_CC
AC_PROG_FC

dnl Find libtool.
AC_PROG_LIBTOOL

dnl Check for assert macros.
AC_HEADER_ASSERT

dnl Check for math library.
AC_CHECK_HEADERS([math.h])
AC_CHECK_LIB([m], [ceil],
             [],
             [AC_MSG_FAILURE([ceil() not found])],
             [])

dnl Check for stdarg.h.
AC_CHECK_HEADERS([stdarg.h])

dnl Write results.
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
Makefile
src/Makefile
tests/Makefile
])

AC_OUTPUT
