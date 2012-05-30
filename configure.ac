AC_INIT([SpAMMPack],
        [beta-2012-05-22],
        [freeon-devel@nongnu.org],
        [],
        [http://www.freeon.org])

AC_COPYRIGHT([
This code is part of the FreeON suite of programs for linear scaling
electronic structure theory and ab initio molecular dynamics.

Copyright (2004). The Regents of the University of California. This
material was produced under U.S. Government contract W-7405-ENG-36
for Los Alamos National Laboratory, which is operated by the University
of California for the U.S. Department of Energy. The U.S. Government has
rights to use, reproduce, and distribute this software.  NEITHER THE
GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 3 of the License, or (at your
option) any later version. Accordingly, this program is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE. See the GNU General Public License at www.gnu.org for details.

While you may do as you like with this software, the GNU license requires
that you clearly mark derivative software.  In addition, you are encouraged
to return derivative works to the FreeON group for review, and possible
disemination in future releases.
])

dnl Initialize automake.
AM_INIT_AUTOMAKE([dist-bzip2 no-dist-gzip silent-rules])
AM_SILENT_RULES([yes])

dnl Get compilers.
AC_PROG_CC
AC_PROG_FC

dnl Check for default compiler flags.
case ${FC} in
  *gfortran*)
    AC_MSG_NOTICE([Detected gfortran, the gcc Fortran compiler])
    if test "${ac_cv_env_FCFLAGS_set}" != "set"; then
      FCFLAGS="-O3 -g -Wall -DTEST_REPEAT=1"
      scaling_reference_file="freeon-10.gcc.reference"
      AC_MSG_NOTICE([setting default compiler flags: ${FCFLAGS}])
    fi
    ;;
  *ifort*)
    AC_MSG_NOTICE([Detected ifort, the Intel Fortran compiler])
    if test "${ac_cv_env_FCFLAGS_set}" != "set"; then
      FCFLAGS="-O3 -no-prec-div -xHost -heap-arrays -g -check all -traceback -DTEST_REPEAT=100"
      scaling_reference_file="freeon-10.Intel.reference"
      AC_MSG_NOTICE([setting default compiler flags: ${FCFLAGS}])
    fi
    ;;
  *) AC_MSG_FAILURE([Undetected FC compiler (${FC})])
     ;;
esac

dnl Set the proper scaling reference file.
AC_SUBST([SCALING_REFERENCE], [${scaling_reference_file}], [The reference file])

dnl Find out how to enable OpenMP in C.
AC_LANG(C)
AC_OPENMP

dnl Find out how to enable OpenMP in Fortran.
AC_LANG(Fortran)
AC_OPENMP

dnl Switch back to C.
AC_LANG(C)

dnl Get libtool.
AC_PROG_LIBTOOL

dnl Define which timer to use.
dnl AC_DEFINE([TIMER_GETRUSAGE], [1], [Use getrusage() to time the library])
AC_DEFINE([TIMER_GETTIMEOFDAY], [1], [Use gettimeofday() to time the library])

dnl Output files.
AC_CONFIG_FILES([
Doxyfile
Makefile
src/Makefile
tests/Makefile
tests/scaling_test.sh
])

AC_OUTPUT
