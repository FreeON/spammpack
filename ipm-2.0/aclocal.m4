dnl ****************************************************************
dnl A combination of AC_CHECK_TOOL and AC_PROG_CC that checks for cross
dnl compilers
dnl ****************************************************************

AC_DEFUN(ACSPERT_TOOL_CC,
[AC_BEFORE([$0], [AC_PROG_CPP])dnl
AC_CHECK_TOOL(CC, gcc, cc)

AC_MSG_CHECKING(whether we are using GNU C)
AC_CACHE_VAL(ac_cv_prog_gcc,
[dnl The semicolon is to pacify NeXT's syntax-checking cpp.
cat > conftest.c <<EOF
#ifdef __GNUC__
  yes;
#endif
EOF
if ${CC-cc} -E conftest.c 2>&AC_FD_CC | egrep yes >/dev/null 2>&1; then
  ac_cv_prog_gcc=yes
else
  ac_cv_prog_gcc=no
fi])dnl
AC_MSG_RESULT($ac_cv_prog_gcc)
if test $ac_cv_prog_gcc = yes; then
  GCC=yes
  if test "${CFLAGS+set}" != set; then
    AC_MSG_CHECKING(whether ${CC-cc} accepts -g)
AC_CACHE_VAL(ac_cv_prog_gcc_g,
[echo 'void f(){}' > conftest.c
if test -z "`${CC-cc} -g -c conftest.c 2>&1`"; then
  ac_cv_prog_gcc_g=yes
else
  ac_cv_prog_gcc_g=no
fi
rm -f conftest*
])dnl
    AC_MSG_RESULT($ac_cv_prog_gcc_g)
    if test $ac_cv_prog_gcc_g = yes; then
      CFLAGS="-g -O"
    else
      CFLAGS="-O"
    fi
  fi
else
  GCC=
  test "${CFLAGS+set}" = set || CFLAGS="-g"
fi
])

dnl ****************************************************************
dnl Work out how to run programs - empty command on a workstation,
dnl but an appropriate spertrun command for SPERT
dnl ****************************************************************
AC_DEFUN(ACSPERT_PROG_RUN,
[
RUN=
AC_MSG_CHECKING(how to run test programs)
if test "$host_cpu" = "torrent0"; then
  RUN=spertrun
fi
AC_SUBST(RUN)
AC_MSG_RESULT(${RUN:-host})
])
