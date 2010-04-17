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

dnl Check for bzip2 library (for tests).
AC_CHECK_HEADERS([bzlib.h])
AC_CHECK_LIB([bz2], [BZ2_bzopen], [], [], [])

dnl Check for blas.
AC_CHECK_LIB([blas], [dgemm],
             [LIBS="-lblas ${LIBS}"
              AC_DEFINE([DGEMM], [dgemm], [We have dgemm().])],
              [], [])
AC_CHECK_LIB([blas], [dgemm_],
             [LIBS="-lblas ${LIBS}"
              AC_DEFINE([DGEMM], [dgemm_], [We have dgemm_().])],
              [], [])
AC_CHECK_LIB([blas], [dgemm__],
             [LIBS="-lblas ${LIBS}"
              AC_DEFINE([DGEMM], [dgemm__], [We have dgemm__().])],
              [], [])

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
