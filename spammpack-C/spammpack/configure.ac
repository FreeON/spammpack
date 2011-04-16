dnl vim:tw=0

dnl Define some basic information of this package.
AC_INIT([spammpack], [2011-02-03], [nicolasbock@gmail.com], [spammpack], [http://freeon.org])

dnl Check for build host.
AC_CANONICAL_TARGET
if test "${target_cpu}" != "x86_64"; then
  AC_MSG_FAILURE([we currently only support the 64-bit x86 architecture])
fi

dnl Initialize automake.
AM_INIT_AUTOMAKE([dist-bzip2 no-dist-gzip silent-rules])
AM_SILENT_RULES([yes])

dnl Define location of m4 macros.
AC_CONFIG_MACRO_DIR([m4])

dnl Identify compilers.
AC_PROG_CC_C99
AC_PROG_F77
AM_PROG_AS
AC_PROG_FC

dnl Find python.
AM_PATH_PYTHON

dnl Find the python include path.
AC_MSG_CHECKING([for python include path])
python_error="no"
python_include=`${PYTHON} -c "from distutils import sysconfig; print sysconfig.get_python_inc()"` || python_error="yes"
if test "${python_error}" = "yes"; then
  AC_MSG_RESULT([an error occured])
  AC_MSG_FAILURE([I need the python include path.])
fi
AC_MSG_RESULT([${python_include}])
AC_SUBST(PYTHON_INCLUDE, ${python_include})

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

dnl Check for SSE.
AC_CHECK_HEADERS([xmmintrin.h])
if test "${ac_cv_header_xmmintrin_h}" = "yes"; then
  AC_DEFINE([HAVE_SSE], [1], [Use SSE instructions in some parts of the code.])
fi

dnl Check for PAPI.
AC_MSG_CHECKING([whether to enable PAPI profiling interface])
AC_ARG_ENABLE([PAPI], [AS_HELP_STRING([--enable-PAPI],
              [enables the PAPI profiling interface])],
              [], [enable_PAPI="no"])

if test "${enable_PAPI}" = "yes"; then
  AC_MSG_RESULT([yes])

  dnl We need to define the C99 macro to get the declaration of caddr_t.
  CPPFLAGS="${CPPFLAGS} -DC99"

  have_PAPI="yes"
  AC_CHECK_HEADERS([papi.h], [], [unset have_PAPI])
  AC_CHECK_LIB([papi], [PAPI_library_init], [], [unset have_PAPI])
  if test "${have_PAPI}" = "yes"; then
    AC_DEFINE([HAVE_PAPI], [1], [Should we use the PAPI profiling interface?])
  else
    AC_MSG_FAILURE([PAPI requested but not available. You might have to
specify the PAPI header file location by adding to CPPFLAGS and/or the PAPI
library location by adding to LDFLAGS.])
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

dnl Define basic tree properties.
have_N_block=4
have_N_child=2
have_N_stride=1
have_kernel_depth=2
have_N_kernel_block=$((${have_N_child}**${have_kernel_depth}))
have_N_kernel=$((${have_N_block}*${have_N_child}**${have_kernel_depth}))

AC_MSG_NOTICE([at the kernel tier, the dense matrix block has a size of ${have_N_kernel}x${have_N_kernel} elements, representing a ${have_N_kernel_block}x${have_N_kernel_block} matrix of ${have_N_block}x${have_N_block} matrix blocks])

dnl Define default memory alignment.
have_alignment=64

dnl Define some constants.
AC_DEFINE_UNQUOTED(SPAMM_ALIGNMENT, [${have_alignment}], [Define the alignment of memory allocations.])
AC_DEFINE_UNQUOTED(SPAMM_PAGE_ALIGNMENT, [4096], [Define the alignment on a 4kB page.])
AC_DEFINE_UNQUOTED(SPAMM_N_BLOCK, [${have_N_block}], [Define the size of a basic matrix block as NxN.])
AC_DEFINE_UNQUOTED(SPAMM_N_CHILD, [${have_N_child}], [Define the size of the children node matrix as NxN.])
AC_DEFINE_UNQUOTED(SPAMM_KERNEL_DEPTH, [${have_kernel_depth}], [Define the number of tiers that are allocated contiguously for the computational kernel.])
AC_DEFINE_UNQUOTED(SPAMM_N_KERNEL_BLOCK, [${have_N_kernel_block}], [Define the size of the kernel matrix as NxN basic kernel blocks.])
AC_DEFINE_UNQUOTED(SPAMM_N_KERNEL, [${have_N_kernel}], [Define the size of the kernel matrix as NxN single elements.])
AC_DEFINE_UNQUOTED(SPAMM_N_STRIDE, [${have_N_stride}], []Define the stride in the kernel.)
AC_DEFINE(SPAMM_BLOCK_ORDER, [1], [Define the block order of the basic matrix blocks in the kernel tier. Possible values are 1 - row major order, 2 - Z-curve order.])

dnl For benchmarking.
AC_DEFINE([SPAMM_MULTIPLY_BETA], [1], [Multiply C matrix with beta factor.])
AC_DEFINE([SPAMM_MULTIPLY_SORT_INDEX], [1], [Sort A and B linear indices.])
AC_DEFINE([SPAMM_MULTIPLY_K_LOOKUP], [1], [Create k lookup table.])
AC_DEFINE([SPAMM_MULTIPLY_SORT_NORM], [1], [Subsort A and B by norms.])
AC_DEFINE([SPAMM_MULTIPLY_COPY_INDICES], [1], [Copy indices to array.])
AC_DEFINE([SPAMM_MULTIPLY_COPY_3D], [1], [Copy 3D indices to array.])
AC_DEFINE([SPAMM_MULTIPLY_CONVOLUTE], [1], [Convolute.])
AC_DEFINE([SPAMM_MULTIPLY_FREE], [1], [Free memory.])
AC_DEFINE([SPAMM_MULTIPLY_STREAM], [1], [Stream multiply.])
AC_DEFINE([SPAMM_MULTIPLY_PRODUCT_COUNT], [1], [Count the number of block products. This is only meant for debugging since it degrades the performance of the multiply.])
AC_DEFINE([SPAMM_MULTIPLY_FINAL_FREE], [1], [Free remaining memory before exit.])

dnl Load git commit ID and set.
if test -f ${srcdir}/localversion; then
  . ${srcdir}/localversion
else
  branch="unknown"
  localversion="unknown"
fi

dnl Fix the package version in the documentation.
AC_SUBST(PACKAGE_VERSION_STRING, [${PACKAGE_VERSION}])
AC_SUBST(LOCAL_VERSION, [${branch}:${localversion}])
AC_DEFINE_UNQUOTED(SPAMM_VERSION, ["${VERSION}"], [The spammpack version.])

dnl Define the git commit tag.
AC_DEFINE_UNQUOTED(SPAMM_COMMIT_TAG, ["${localversion}"], [Define the git commit tag of this build.])

dnl Some output.
AC_MSG_NOTICE([configuring local version ${branch}:${localversion}])

dnl Store the command line.
AC_DEFINE_UNQUOTED([SPAMM_CONFIGURE_ARGUMENTS], ["${ac_configure_args}"], [The command line with which spammpack was configured.])

dnl Write results.
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
Doxyfile
Makefile
interfaces/Makefile
interfaces/f90/Makefile
interfaces/python/Makefile
src/Makefile
src/spamm_manual_mainpage.h
tests/Makefile
])

AC_CONFIG_COMMANDS([src/spamm_config.h], [
  echo -n "configure: creating src/spamm_config.h... "
  echo "/** @file */" > src/spamm_config.h.temp
  echo >> src/spamm_config.h.temp
  echo "#ifndef __SPAMM_CONFIG_H" >> src/spamm_config.h.temp
  echo "#define __SPAMM_CONFIG_H" >> src/spamm_config.h.temp
  echo >> src/spamm_config.h.temp
  echo "/* This file contains the configuration options that directly" >> src/spamm_config.h.temp
  echo " * affect the built of spamm." >> src/spamm_config.h.temp
  echo " */" >> src/spamm_config.h.temp
  echo >> src/spamm_config.h.temp
  grep "#define SPAMM" config.h >> src/spamm_config.h.temp
  echo >> src/spamm_config.h.temp
  echo "#endif" >> src/spamm_config.h.temp
  if (diff --brief src/spamm_config.h src/spamm_config.h.temp > /dev/null 2>&1); then
    echo "no difference, keeping old one"
    rm -f src/spamm_config.h.temp
  else
    echo "updating src/spamm_config.h"
    mv -f src/spamm_config.h.temp src/spamm_config.h
  fi
])

AC_OUTPUT
