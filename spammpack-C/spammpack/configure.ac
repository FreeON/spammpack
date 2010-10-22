AC_INIT([spammpack], [2010-06-24], [nicolasbock@gmail.com])
AM_INIT_AUTOMAKE

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

dnl If `stdint.h' or `inttypes.h' defines the type `intptr_t', define
dnl `HAVE_INTPTR_T'.  Otherwise, define `intptr_t' to a signed integer type wide
dnl enough to hold a pointer, if such a type exists.
AC_TYPE_INTPTR_T

dnl Define `pid_t' to a suitable type, if standard headers do not define it.
AC_TYPE_PID_T

dnl Define precision.
AC_ARG_WITH([precision],
            [AS_HELP_STRING([--with-precision], [set the precision of floating
             point numbers to single or double [default=single]])],
             [], [with_precision="single"])

if test "${with_precision}" = "single"; then
  AC_MSG_NOTICE([using single precision])
  AC_DEFINE([FLOATING_PRECISION], [float], [Define the floating point number type])
elif test "${with_precision}" = "double"; then
  AC_MSG_NOTICE([using double precision])
  AC_DEFINE([FLOATING_PRECISION], [double], [Define the floating point number type])
else
  AC_MSG_FAILURE([unknown precision: ${with_precision}])
fi

dnl Check for math library.
AC_CHECK_HEADERS([math.h])
AC_CHECK_LIB([m], [ceil],
             [],
             [AC_MSG_FAILURE([ceil() not found])],
             [])

dnl Check for malloc with alignment.
AC_DEFINE([_POSIX_C_SOURCE], [200809L], [Define POSIX to enable
           posix_memalign()])
AC_CHECK_DECL([posix_memalign],
              [declared_posix_memalign="yes"],
              [declared_posix_memalign="no"], [])

if test "${declared_posix_memalign}" = "yes"; then
  AC_CHECK_LIB([c], [posix_memalign],
               [AC_DEFINE([HAVE_POSIX_MEMALIGN], [1], [We can allocate aligned
                memory])], [], [])
fi

dnl Check for blas.
AC_MSG_CHECKING([whether to use BLAS])
AC_ARG_ENABLE([blas], [AS_HELP_STRING([--enable-blas], [enables the use of BLAS routines])],
              [], [enable_blas="no"])

if test "${enable_blas}" = "yes"; then
  AC_MSG_RESULT([yes])

  if test "${with_precision}" = "double"; then
    AC_CHECK_LIB([blas], [dgemm],
                 [LIBS="-lblas ${LIBS}"
                  have_blas="yes"
                  AC_DEFINE([DGEMM], [dgemm], [How many underscores are in dgemm()?])
                  AC_DEFINE([HAVE_BLAS], [1], [whether we can use blas routines])],
                  [], [])
    AC_CHECK_LIB([blas], [dgemm_],
                 [LIBS="-lblas ${LIBS}"
                  have_blas="yes"
                  AC_DEFINE([DGEMM], [dgemm_], [How many underscores are in dgemm()?])
                  AC_DEFINE([HAVE_BLAS], [1], [whether we can use blas routines])],
                  [], [])
    AC_CHECK_LIB([blas], [dgemm__],
                 [LIBS="-lblas ${LIBS}"
                  have_blas="yes"
                  AC_DEFINE([DGEMM], [dgemm__], [How many underscores are in dgemm()?])
                  AC_DEFINE([HAVE_BLAS], [1], [whether we can use blas routines])],
                  [], [])
  elif test "${with_precision}" = "single"; then
    AC_CHECK_LIB([blas], [sgemm],
                 [LIBS="-lblas ${LIBS}"
                  have_blas="yes"
                  AC_DEFINE([DGEMM], [sgemm], [How many underscores are in sgemm()?])
                  AC_DEFINE([HAVE_BLAS], [1], [whether we can use blas routines])],
                  [], [])
    AC_CHECK_LIB([blas], [sgemm_],
                 [LIBS="-lblas ${LIBS}"
                  have_blas="yes"
                  AC_DEFINE([DGEMM], [sgemm_], [How many underscores are in sgemm()?])
                  AC_DEFINE([HAVE_BLAS], [1], [whether we can use blas routines])],
                  [], [])
    AC_CHECK_LIB([blas], [sgemm__],
                 [LIBS="-lblas ${LIBS}"
                  have_blas="yes"
                  AC_DEFINE([DGEMM], [sgemm__], [How many underscores are in sgemm()?])
                  AC_DEFINE([HAVE_BLAS], [1], [whether we can use blas routines])],
                  [], [])
  else
    AC_MSG_FAILURE([unknown precision])
  fi
else
  AC_MSG_RESULT([no])
fi

dnl Check for whether to use dummy kernel. This can be useful to time the tree
dnl work overhead.
AC_MSG_CHECKING([whether to use dummy kernel])
AC_ARG_ENABLE([dummy-kernel], [AS_HELP_STRING([--enable-dummy-kernel],
              [enables the dummy matrix-matrix multiply kernel])],
              [], [enable_dummy_kernel="no"])

if test "${enable_dummy_kernel}" = "yes"; then
  AC_MSG_RESULT([yes])
  have_blas="yes"
  AC_DEFINE([DGEMM], [spamm_sgemm_dummy], [The dummy matrix kernel])
else
  AC_MSG_RESULT([no])
fi

AM_CONDITIONAL([HAVE_BLAS], [test "${have_blas}" = "yes"])

dnl Check for location of cuda.
AC_MSG_CHECKING([whether to enable cuda])
AC_ARG_ENABLE([cuda], [AS_HELP_STRING([--enable-cuda=path_to_sdk], [Enable
               cuda and set the location of the cuda sdk. /include and /lib
               will be appended to path_to_sdk to find headers and
               libraries.])],
               [], [enable_cuda="no"])

if test "${enable_cuda}" != "no"; then
  AC_MSG_RESULT([yes, using cuda path ${enable_cuda}])
  have_cuda="yes"
  old_CPPFLAGS=${CPPFLAGS}
  CPPFLAGS="${CPPFLAGS} -I${enable_cuda}/include"
  AC_CHECK_HEADERS([cublas.h], [],
                   [CPPFLAGS=${old_CPPFLAGS}
                    unset have_cuda])
  old_LDFLAGS=${LDFLAGS}
  LDFLAGS="${LDFLAGS} -L${enable_cuda}/lib64 -L${enable_cuda}/lib"
  #LDFLAGS="${LDFLAGS} -L${enable_cuda}/lib64"
  AC_CHECK_LIB([cublas], [cublasSgemm], [],
               [LDFLAGS=${old_LDFLAGS}
                unset have_cuda], [])
  unset old_CPPFLAGS
  unset old_LDFLAGS
  if test "${have_cuda}" = "yes"; then
    AC_DEFINE([HAVE_CUDA], [1], [whether we can use cuda])
  fi
else
  AC_MSG_RESULT([no])
fi

AM_CONDITIONAL([HAVE_CUDA], [test "${have_cuda}" = "yes"])

dnl Check for stdarg.h.
AC_CHECK_HEADERS([stdarg.h])

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
    AC_DEFINE([HAVE_PAPI], [1], [whether we should use the PAPI profiling interface])
  else
    AC_MSG_FAILURE([PAPI requested but not available. You might have to
                    specify the PAPI header file location by adding to
                    CPPFLAGS the PAPI library location by adding to LDFLAGS.])
  fi
else
  AC_MSG_RESULT([no])
fi

dnl Check for SSE intrinsics.
AC_CHECK_HEADERS([xmmintrin.h])

dnl Check for pthread support.
have_PTHREAD="yes"
AC_CHECK_HEADERS([pthread.h], [], [unset have_PTHREAD])
AC_CHECK_LIB([pthread], [pthread_create], [], [unset have_PTHREAD])
if test "${have_PTHREAD}" = "yes"; then
  AC_DEFINE([HAVE_PTHREAD], [1], [whether we have pthread support])
else
  AC_MSG_FAILURE([no pthread support])
fi

dnl Check for memory alignment.
AC_MSG_CHECKING([for memory alignment])
AC_ARG_WITH([alignment], [AS_HELP_STRING([--with-alignment], [Define the
             alignment of memory allocations. The default is 64, a cache
             line.])],
             [have_alignment="${withval}"],
             [have_alignment="64"])
AC_MSG_RESULT([alignment on ${have_alignment} byte boundaries])

dnl Check for matrix block size.
AC_MSG_CHECKING([for number of rows of matrix block])
AC_ARG_WITH([N-block], [AS_HELP_STRING([--with-N-block], [Define the size of a
             matrix block as NxN. The default is N = 4, i.e. matrix blocks of
             size 4x4.])],
             [have_N_block="${withval}"],
             [have_N_block="4"])
AC_MSG_RESULT([using ${have_N_block}x${have_N_block} matrix blocks])

dnl Check for children node matrix.
AC_MSG_CHECKING([for number of columns of children node matrix])
AC_ARG_WITH([N-child], [AS_HELP_STRING([--with-N-child], [Define the size of
             the children node matrix as NxN. The default is N = 2, i.e. a 2x2
             children node matrix, which defines a quadtree.])],
             [have_N_child="${withval}"],
             [have_N_child="2"])
AC_MSG_RESULT([using a ${have_N_child}x${have_N_child} children node matrix])

dnl Check for kernel depth.
AC_MSG_CHECKING([for depth in tiers of computational kernel])
AC_ARG_WITH([kernel-depth], [AS_HELP_STRING([--with-kernel-depth], [Define
             the kernel depth in terms of tiers that allocate the matrix
             blocks contiguously for the computational kernel. This value has
             to match the computational kernel.])],
             [have_kernel_depth="${withval}"],
             [have_kernel_depth="2"])
AC_MSG_RESULT([using ${have_kernel_depth} tier(s)])

dnl Define some constants.
AC_DEFINE_UNQUOTED([SPAMM_ALIGNMENT], [${have_alignment}], [Define the alignment of memory allocations])
AC_DEFINE_UNQUOTED([SPAMM_N_BLOCK], [${have_N_block}], [Define the size of a matrix block as NxN])
AC_DEFINE_UNQUOTED([SPAMM_N_CHILD], [${have_N_child}], [Define the size of the children node matrix as NxN])
AC_DEFINE_UNQUOTED([SPAMM_KERNEL_DEPTH], [${have_kernel_depth}], [Define the number of tiers that are allocated contiguously for the computational kernel])
AC_DEFINE_UNQUOTED([SPAMM_N_KERNEL], [$((${have_N_block}*${have_N_child}**${have_kernel_depth}))], [Define the size of the matrix kernel as NxN])

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
AC_MSG_NOTICE([configuring local version ${branch}:${localversion}])

dnl Write results.
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
Doxyfile
Makefile
benchmarks/Makefile
figures/Makefile
prototype/Makefile
prototype/Z_curve/Makefile
prototype/cublas/Makefile
prototype/cuda_gemm_kernel/Makefile
prototype/kernel_comparison/Makefile
prototype/linear_quadtree/Makefile
prototype/stream_multiply/Makefile
sparsekit/Makefile
src/Makefile
src/spamm_mainpage.h
tests/Makefile
tests/add/Makefile
tests/ll/Makefile
tests/mm/Makefile
tests/multiply/Makefile
])

AC_OUTPUT
