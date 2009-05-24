AC_INIT([LinearAlgebraLibrary], [0.1], [nbock@lanl.gov])
AM_INIT_AUTOMAKE

dnl Identify compilers.
AC_PROG_CC
AC_PROG_FC

dnl Find libtool.
AC_PROG_LIBTOOL

dnl Check for assert macros.
AC_HEADER_ASSERT

dnl Get pointer length.
AC_LANG(C)
AC_MSG_CHECKING([for pointer length])

pointer_length=-1

if test ${pointer_length} -lt 0; then

  AC_RUN_IFELSE([AC_LANG_SOURCE([int main ()
{
  if (sizeof(void *) == 4) { return 0; }
  else { return -1; }
}])],
                [pointer_length=4], [])

fi

if test ${pointer_length} -lt 0; then

  AC_RUN_IFELSE([AC_LANG_SOURCE([int main ()
{
  if (sizeof(void *) == 8) { return 0; }
  else { return -1; }
}])],
                [pointer_length=8], [])

fi

if test ${pointer_length} -lt 0; then
  AC_MSG_RESULT([unknown])
  AC_MSG_FAILURE([can not determine pointer length])
fi

AC_MSG_RESULT(${pointer_length})
AC_DEFINE_UNQUOTED([POINTER_LENGTH], [${pointer_length}], [The length of a pointer])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
Makefile
interface/Makefile
interface/f90/Makefile
src/Makefile
tests/Makefile
])

AC_OUTPUT
