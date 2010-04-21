AC_INIT([SpAMM], [2010-04-19], [nbock@lanl.gov])
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

dnl Check for blas.
AC_MSG_CHECKING([whether to use BLAS])
AC_ARG_ENABLE([blas], [AS_HELP_STRING([--enable-blas], [enables the use of BLAS routines])],
              [], [enable_blas="no"])

if test "${enable_blas}" = "yes"; then
  AC_MSG_RESULT([yes])
  AC_CHECK_LIB([blas], [dgemm],
               [LIBS="-lblas ${LIBS}"
                AC_DEFINE([DGEMM], [dgemm], [How many underscores are in dgemm()?])],
                [], [])
  AC_CHECK_LIB([blas], [dgemm_],
               [LIBS="-lblas ${LIBS}"
                AC_DEFINE([DGEMM], [dgemm_], [How many underscores are in dgemm()?])],
                [], [])
  AC_CHECK_LIB([blas], [dgemm__],
               [LIBS="-lblas ${LIBS}"
                AC_DEFINE([DGEMM], [dgemm__], [How many underscores are in dgemm()?])],
                [], [])
else
  AC_MSG_RESULT([no])
fi

dnl Check for stdarg.h.
AC_CHECK_HEADERS([stdarg.h])

dnl Write results.
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
Makefile
src/Makefile
tests/Makefile
benchmarks/Makefile
prototype/Makefile
prototype/cublas/Makefile
])

AC_OUTPUT
