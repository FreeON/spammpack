AC_INIT([migration], [1], [nicolasbock@freeon.org])
AM_INIT_AUTOMAKE([dist-bzip2 no-dist-gzip silent-rules])
AM_SILENT_RULES([yes])

AC_ARG_VAR([CHARMC], [Charm++ compiler command])
if test -z "${CHARMC}"; then
  CHARMC="charmc"
fi
CXX="${CHARMC}"

AC_PROG_CC
AC_PROG_CXX

dnl Test everything else with the C compiler. The CHARMC wrapper can not be
dnl tested with standard autoconf macros.
AC_LANG(C)

dnl Test for OpenMP flags.
AC_OPENMP
CXXFLAGS="${CXXFLAGS:+${CXXFLAGS} }${OPENMP_CFLAGS}"

dnl Test for some libraries.
AC_SEARCH_LIBS([sqrt], [m])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
