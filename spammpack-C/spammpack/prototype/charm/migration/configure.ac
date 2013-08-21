AC_INIT([migration], [1], [nicolasbock@freeon.org])
AM_INIT_AUTOMAKE([dist-bzip2 no-dist-gzip silent-rules])
AM_SILENT_RULES([yes])

AC_ARG_VAR([CHARMC], [Charm++ compiler command])
if test -z "${CHARMC}"; then
  CHARMC="charmc"
fi
AC_CHECK_PROG([CHARMC],
              ["${CHARMC}"],
              ["${CHARMC}"],
              [undefined],
              [PATH = $(dirname "${CHARMC}") ])
if test "${CHARMC}" = "undefined"; then
  AC_MSG_FAILURE([need CHARMC])
fi

AC_PROG_CC
AC_PROG_CXX([${CHARMC}])

AC_LANG(C++)

AC_SEARCH_LIBS([sqrt], [m])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
