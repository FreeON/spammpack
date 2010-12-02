dnl vim:tw=0

dnl Define some basic information of this package.
AC_INIT([spammpack], [2010-11-29], [nicolasbock@gmail.com], [spammpack], [http://freeon.org])
AM_INIT_AUTOMAKE([dist-bzip2 no-dist-gzip silent-rules])
AM_SILENT_RULES([yes])

dnl Define location of m4 macros.
AC_CONFIG_MACRO_DIR([m4])

dnl Identify compilers.
AC_PROG_CC_STDC
AC_PROG_F77
AM_PROG_AS

dnl Find libtool.
AC_PROG_LIBTOOL

dnl Check for assert macros.
AC_HEADER_ASSERT

dnl Check for glib.
glib_libs=( "glib-2.0" )
glib_libs_ac=( "glib_2_0" )
glib_cppflags=( "" "-I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include" )
old_CPPFAGS=$CPPFLAGS
old_LIBS=$LIBS

for (( i = 0; i < ${#glib_libs[*]}; i++ )); do
  for (( j = 0; j < ${#glib_cppflags[*]}; j++ )); do

    dnl We have to use Quadrigraphs to write '[' and ']'. Otherwise m4 will not
    dnl print the brackets. The quadrigraphs for the brackets are:
    dnl
    dnl "@<:@" = "[" and "@:>@" = "]".
    have_GLIB="yes"
    unset ac_cv_header_glib_h
    unset ac_cv_lib_${glib_libs_ac@<:@$i@:>@}___g_hash_table_new
    CPPFLAGS="${old_CPPFLAGS} ${glib_cppflags@<:@$j@:>@}"
    LIBS=$old_LIBS
    AC_MSG_NOTICE([checking for ${glib_libs@<:@$i@:>@}, using CPPFLAGS="${CPPFLAGS}"])
    AC_CHECK_HEADERS([glib.h], [], [unset have_GLIB])
    AC_CHECK_LIB([${glib_libs@<:@$i@:>@}], [g_hash_table_new], [], [unset have_GLIB])

    if test "${have_GLIB}" = "yes"; then
      break
    else
      AC_MSG_NOTICE([test for ${glib_libs@<:@$i@:>@} failed])
    fi
  done

  if test "${have_GLIB}" = "yes"; then
    break
  fi
done

if test "${have_GLIB}" != "yes"; then
  AC_MSG_FAILURE([can not find glib. Consider add custom CPPFLAGS, LDFLAGS, and/or LIBS to configure])
else
  AC_MSG_NOTICE([test for ${glib_libs@<:@$i@:>@} passed, found working glib])
fi

dnl Check for math library.
AC_CHECK_HEADERS([math.h])
AC_CHECK_LIB([m], [ceil],
             [],
             [AC_MSG_FAILURE([ceil() not found])],
             [])

dnl Check for malloc with alignment.
AC_DEFINE([_POSIX_C_SOURCE], [200809L], [Define POSIX to enable posix_memalign().])
AC_CHECK_DECL([posix_memalign],
              [declared_posix_memalign="yes"],
              [declared_posix_memalign="no"], [])

if test "${declared_posix_memalign}" = "yes"; then
  AC_CHECK_LIB([c], [posix_memalign],
               [AC_DEFINE([HAVE_POSIX_MEMALIGN], [1], [We can allocate aligned
                memory.])], [], [])
fi

dnl Check for PAPI.
AC_MSG_CHECKING([whether to enable PAPI profiling interface])
AC_ARG_ENABLE([PAPI], [AS_HELP_STRING([--enable-PAPI],
              [enables the PAPI profiling interface])],
              [], [enable_PAPI="no"])

if test "${enable_PAPI}" = "yes"; then
  AC_MSG_RESULT([yes])
  have_PAPI="yes"
  AC_CHECK_HEADERS([papi.h], [], [unset have_PAPI])
  AC_CHECK_LIB([papi], [PAPI_library_init], [], [unset have_PAPI])
  if test "${have_PAPI}" = "yes"; then
    AC_DEFINE([HAVE_PAPI], [1], [Should we use the PAPI profiling interface?])
  else
    AC_MSG_FAILURE([PAPI requested but not available. You might have to
                    specify the PAPI header file location by adding to
                    CPPFLAGS the PAPI library location by adding to LDFLAGS.])
  fi
else
  AC_MSG_RESULT([no])
fi

dnl Check for pthread support.
have_PTHREAD="yes"
AC_CHECK_HEADERS([pthread.h], [], [unset have_PTHREAD])
AC_CHECK_LIB([pthread], [pthread_create], [], [unset have_PTHREAD])
if test "${have_PTHREAD}" = "yes"; then
  AC_DEFINE([HAVE_PTHREAD], [1], [Do we have pthread support?])
else
  AC_MSG_FAILURE([no pthread support])
fi

dnl Check for logging.
AC_MSG_CHECKING([whether to enable logging])
AC_ARG_ENABLE([logging], [AS_HELP_STRING([--disable-logging],
              [disables logging output])],
              [], [enable_logging="yes"])

dnl Check for memory alignment.
AC_MSG_CHECKING([for memory alignment])
AC_ARG_WITH([alignment], [AS_HELP_STRING([--with-alignment], [Define the
             alignment of memory allocations. The default is 64, a cache
             line.])],
             [have_alignment="${withval}"],
             [have_alignment="64"])
AC_MSG_RESULT([alignment on ${have_alignment} byte boundaries])

dnl Check for matrix block size.
AC_MSG_CHECKING([for size of matrix block])
AC_ARG_WITH([N-block], [AS_HELP_STRING([--with-N-block], [Define the size of a
             matrix block as NxN. The default is N = 4, i.e. matrix blocks of
             size 4x4.])],
             [have_N_block="${withval}"],
             [have_N_block="4"])
AC_MSG_RESULT([using ${have_N_block}x${have_N_block} matrix blocks])

dnl Check for children node matrix.
AC_MSG_CHECKING([for size of children node matrix])
AC_ARG_WITH([N-child], [AS_HELP_STRING([--with-N-child], [Define the size of
             the children node matrix as NxN. The default is N = 2, i.e. a 2x2
             children node matrix, which defines a quadtree. You better not
             change this unless you really know what you are doing.])],
             [have_N_child="${withval}"],
             [have_N_child="2"])
AC_MSG_RESULT([using a ${have_N_child}x${have_N_child} children node matrix])

if test ${have_N_child} -eq 1; then
  AC_MSG_FAILURE([It is really pointless to define a tree in which each node has only 1 child.])
fi

if test ${have_N_child} -ne 2; then
  AC_MSG_NOTICE([You requested to not generate a quadtree. Please fasten your seatbelt and hope for a miracle...])
fi

dnl Check for kernel depth.
AC_MSG_CHECKING([for depth in tiers of computational kernel])
AC_ARG_WITH([kernel-depth], [AS_HELP_STRING([--with-kernel-depth], [Define
             the kernel depth in terms of tiers that allocate the matrix
             blocks contiguously for the computational kernel. This value has
             to match the computational kernel.])],
             [have_kernel_depth="${withval}"],
             [have_kernel_depth="2"])
AC_MSG_RESULT([using ${have_kernel_depth} tier(s)])

have_N_kernel_block=$((${have_N_child}**${have_kernel_depth}))
have_N_kernel=$((${have_N_block}*${have_N_child}**${have_kernel_depth}))

AC_MSG_NOTICE([at the kernel tier, the dense matrix block has a size of ${have_N_kernel}x${have_N_kernel} elements, representing a ${have_N_kernel_block}x${have_N_kernel_block} matrix of ${have_N_block}x${have_N_block} matrix blocks])

dnl Define some constants.
AC_SUBST(SPAMM_ALIGNMENT, [${have_alignment}])
AC_SUBST(SPAMM_N_BLOCK, [${have_N_block}])
AC_SUBST(SPAMM_N_CHILD, [${have_N_child}])
AC_SUBST(SPAMM_KERNEL_DEPTH, [${have_kernel_depth}])
AC_SUBST(SPAMM_N_KERNEL_BLOCK, [${have_N_kernel_block}])
AC_SUBST(SPAMM_N_KERNEL, [${have_N_kernel}])

dnl Fix the package version everywhere.
AC_SUBST(PACKAGE_VERSION_STRING, [${PACKAGE_VERSION}])

dnl Load git commit ID and set.
if test -f ${srcdir}/localversion; then
  . ${srcdir}/localversion
else
  branch="unknown"
  localversion="unknown"
fi
AC_SUBST(LOCAL_VERSION, [${branch}:${localversion}])
AC_SUBST(COMMIT_TAG, [${localversion}])
AC_MSG_NOTICE([configuring local version ${branch}:${localversion}])

dnl Write results.
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
Doxyfile
Makefile
src/Makefile
src/spamm_config.h
src/spamm_manual_mainpage.h
tests/Makefile
])

AC_OUTPUT
