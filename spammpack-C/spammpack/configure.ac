dnl vim:tw=0

dnl Define some basic information of this package.
AC_INIT([spammpack],
        [2011-02-03],
        [freeon-devel@nongnu.org],
        [],
        [http://www.freeon.org])

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
AC_PROG_CC_STDC
AC_PROG_F77
AM_PROG_AS
AC_PROG_FC

dnl Find python.
AM_PATH_PYTHON

dnl Find the python include path.
AC_MSG_CHECKING([for python include path])
python_error="no"
python_include=`${PYTHON} -c "from distutils import sysconfig; print(sysconfig.get_python_inc())"` || python_error="yes"
if test "${python_error}" = "yes"; then
  AC_MSG_RESULT([an error occured])
  AC_MSG_FAILURE([I need the python include path.])
fi
AC_MSG_RESULT([${python_include}])
AC_SUBST(PYTHON_INCLUDE, ${python_include})

dnl Find OpenMP flags for C.
AC_LANG_PUSH(C)
AC_OPENMP
AC_LANG_POP(C)

dnl Find OpenMP flags for Fortran.
AC_LANG_PUSH(Fortran)
AC_OPENMP
AC_LANG_POP(Fortran)

dnl Find libtool.
AC_PROG_LIBTOOL

dnl Check for assert macros.
AC_HEADER_ASSERT

dnl Check for math library.
AC_SEARCH_LIBS([ceil], [m], [], [AC_MSG_FAILURE([could not find ceil()])])

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
              [enable_PAPI=$enableval],
              [enable_PAPI="no"])
AC_MSG_RESULT([$enable_PAPI])

if test "${enable_PAPI}" = "yes"; then
  dnl We need to define the C99 macro to get the declaration of caddr_t.
  CPPFLAGS="${CPPFLAGS} -DC99"

  have_PAPI="yes"
  AC_CHECK_HEADERS([papi.h], [], [unset have_PAPI])
  AC_SEARCH_LIBS([PAPI_library_init], [papi], [], [unset have_PAPI])
  if test "${have_PAPI}" = "yes"; then
    AC_DEFINE([HAVE_PAPI], [1], [Should we use the PAPI profiling interface?])
  else
    AC_MSG_FAILURE([PAPI requested but not available. You might have to
specify the PAPI header file location by adding to CPPFLAGS and/or the PAPI
library location by adding to LDFLAGS.])
  fi
fi

dnl Check for pthread support.
have_PTHREAD="yes"
AC_CHECK_HEADERS([pthread.h], [], [unset have_PTHREAD])
AC_SEARCH_LIBS([pthread_create], [pthread], [], [unset have_PTHREAD])
if test "${have_PTHREAD}" = "yes"; then
  AC_DEFINE([HAVE_PTHREAD], [1], [Do we have pthread support?])
else
  AC_MSG_FAILURE([no pthread support])
fi

dnl Check for sgemm_().
AC_SEARCH_LIBS([sgemm_], [],
               [have_sgemm="yes"],
               [have_sgemm="no"])

if test "${have_sgemm}" = "no"; then
  AC_MSG_NOTICE([no external sgemm_(), will use internal sgemm()])
  AC_DEFINE([SGEMM], [spamm_sgemm], [SpAMM internal sgemm.])
fi

dnl Check for internal lapack (in case this package is compiled as part of FreeON).
AC_MSG_CHECKING([for internal lapack])
AC_ARG_ENABLE([internal-lapack],
              [AS_HELP_STRING([--enable-internal-lapack],
                              [enables internal lapack (as part of FreeON)])],
              [enable_internal_lapack=$enableval],
              [enable_internal_lapack="no"])
AC_MSG_RESULT([$enable_internal_lapack])
AM_CONDITIONAL(INTERNAL_LAPACK, [test "${enable_internal_lapack}" = "yes"])

dnl Check for logging.
AC_MSG_CHECKING([whether to enable logging])
AC_ARG_ENABLE([logging], [AS_HELP_STRING([--disable-logging],
              [disables logging output])],
              [enable_logging=$enableval],
              [enable_logging="yes"])
AC_MSG_RESULT([$enable_logging])

dnl Find out whether we support SSE and SSE4.1.
AC_MSG_CHECKING([for SSE support])
AC_RUN_IFELSE([AC_LANG_PROGRAM([],
              [asm("mulps %xmm0, %xmm0");])],
              [have_SSE_in_CPU="yes"],
              [have_SSE_in_CPU="no"])
if test "${have_SSE_in_CPU}" = "yes"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi
AM_CONDITIONAL(HAVE_SSE_IN_CPU, [test "${have_SSE_in_CPU}" = "yes"])

AC_MSG_CHECKING([for SSE4.1 support])
AC_RUN_IFELSE([AC_LANG_PROGRAM([],
              [asm("dpps \$0x0, %xmm0, %xmm0");])],
              [have_SSE4_1_in_CPU="yes"],
              [have_SSE4_1_in_CPU="no"])
if test "${have_SSE4_1_in_CPU}" = "yes"; then
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi
AM_CONDITIONAL(HAVE_SSE4_1_IN_CPU, [test "${have_SSE4_1_in_CPU}" = "yes"])

dnl Define basic tree properties.
have_N_block=4
have_N_stride=1
have_kernel_depth=2
have_N_kernel_blocked=$((2**${have_kernel_depth}))
have_N_kernel=$((${have_N_block}*2**${have_kernel_depth}))

AC_MSG_NOTICE([at the kernel tier, the dense matrix block has a size of ${have_N_kernel}x${have_N_kernel} elements, representing a ${have_N_kernel_block}x${have_N_kernel_block} matrix of ${have_N_block}x${have_N_block} matrix blocks])

dnl Define default memory alignment.
have_alignment=64

dnl Define some constants.
AC_DEFINE_UNQUOTED(SPAMM_ALIGNMENT, [${have_alignment}], [Define the alignment of memory allocations.])
AC_DEFINE_UNQUOTED(SPAMM_PAGE_ALIGNMENT, [4096], [Define the alignment on a 4kB page.])
AC_DEFINE_UNQUOTED(SPAMM_N_BLOCK, [${have_N_block}], [Define the size of a basic matrix block as NxN single matrix elements.])
AC_DEFINE_UNQUOTED(SPAMM_KERNEL_DEPTH, [${have_kernel_depth}], [Define the number of tiers that are allocated contiguously for the computational kernel.])
AC_DEFINE_UNQUOTED(SPAMM_N_KERNEL_BLOCKED, [${have_N_kernel_blocked}], [Define the size of the kernel matrix as NxN basic kernel blocks.])
AC_DEFINE_UNQUOTED(SPAMM_N_KERNEL, [${have_N_kernel}], [Define the size of the kernel matrix as NxN single elements.])
AC_DEFINE_UNQUOTED(SPAMM_N_STRIDE, [${have_N_stride}], []Define the stride in the kernel.)
AC_DEFINE(SPAMM_BLOCK_ORDER, [1], [Define the block order of the basic matrix blocks in the kernel tier. Possible values are 1 - row major order, 2 - Z-curve order.])

dnl For benchmarking.
#AC_DEFINE([SPAMM_MULTIPLY_BETA], [1], [Multiply C matrix with beta factor.])
AC_DEFINE([SPAMM_MULTIPLY_SORT_INDEX], [1], [Sort A and B linear indices.])
AC_DEFINE([SPAMM_MULTIPLY_K_LOOKUP], [1], [Create k lookup table.])
AC_DEFINE([SPAMM_MULTIPLY_SORT_NORM], [1], [Subsort A and B by norms.])
AC_DEFINE([SPAMM_MULTIPLY_COPY_INDICES], [1], [Copy indices to array.])
AC_DEFINE([SPAMM_MULTIPLY_CONVOLUTE], [1], [Convolute.])
AC_DEFINE([SPAMM_MULTIPLY_FREE], [1], [Free memory.])
AC_DEFINE([SPAMM_MULTIPLY_STREAM], [1], [Stream multiply.])
AC_DEFINE([SPAMM_MULTIPLY_FINAL_FREE], [1], [Free remaining memory before exit.])
AC_DEFINE([SPAMM_MULTIPLY_UPDATE_NORM], [1], [Update matrix norms on C tree after multiply.])
#AC_DEFINE([SPAMM_MULTIPLY_PRINT_ALOT], [1], [Print a lot of information while multiplying.])

AC_MSG_CHECKING([whether to count the number of basic matrix block products])
AC_ARG_ENABLE([product-count],
              [AS_HELP_STRING([--enable-product-count], [Enable basic matrix block product count.])],
              [use_product_count=${enableval}],
              [use_product_count="no"])
AC_MSG_RESULT([$use_product_count])

if test "${use_product_count}" = "yes"; then
  AC_DEFINE([SPAMM_MULTIPLY_PRODUCT_COUNT], [1], [Count the number of block products. This is only meant for debugging since it degrades the performance of the multiply.])
fi

AC_MSG_CHECKING([whether to use real hashtables])
AC_ARG_ENABLE([hashtable],
              [AS_HELP_STRING([--enable-hashtable], [Enables the use of real hashtables as opposed to arrays.])],
              [use_hashtable=${enableval}],
              [use_hashtable="no"])
AC_MSG_RESULT([$use_hashtable])

use_hashtable_hash=0
if test "${use_hashtable}" = "yes"; then
  dnl Use hashtables in matrix tiers?
  AC_DEFINE([SPAMM_USE_HASHTABLE], [1], [Use tier hashtables.])

  AC_MSG_CHECKING([whether to use the direct hash])
  AC_ARG_ENABLE([hash-direct],
                [AS_HELP_STRING([--enable-hash-direct], [Enables the direct hash, i.e. the key is not modified.])],
                [use_hashtable_hash_direct=${enableval}],
                [use_hashtable_hash_direct="no"])
  AC_MSG_RESULT([$use_hashtable_hash_direct])

  AC_MSG_CHECKING([whether to use the Jenkins hash])
  AC_ARG_ENABLE([hash-jenkins],
                [AS_HELP_STRING([--enable-hash-jenkins], [Enables the lookup3 hash by Jenkins])],
                [use_hashtable_hash_jenkins=${enableval}],
                [use_hashtable_hash_jenkins="no"])
  AC_MSG_RESULT([$use_hashtable_hash_jenkins])

  AC_MSG_CHECKING([whether to use the MurmurHash3 hash])
  AC_ARG_ENABLE([hash-murmurhash3],
                [AS_HELP_STRING([--enable-hash-murmurhash3], [Enables the MurmurHash3 hash])],
                [use_hashtable_hash_murmurhash3=${enableval}],
                [use_hashtable_hash_murmurhash3="no"])
  AC_MSG_RESULT([$use_hashtable_hash_murmurhash3])

  dnl Set default.
  use_hashtable_hash=2

  dnl Define hashtable name.
  current_hashtable_name="Jenkins"

  dnl Test what the user wanted.
  if test "$use_hashtable_hash_direct" = "yes"; then
    use_hashtable_hash=1
  elif test "$use_hashtable_hash_jenkins" = "yes"; then
    use_hashtable_hash=2
  elif test "$use_hashtable_hash_murmurhash3" = "yes"; then
    use_hashtable_hash=3
  else
    AC_MSG_NOTICE([no hash preference given, using $current_hashtable_name])
  fi
fi
AC_DEFINE_UNQUOTED([SPAMM_HASHTABLE_HASH], [$use_hashtable_hash], [The type of hash function])

AC_MSG_CHECKING([what convolution implementation to use])
AC_ARG_WITH([convolution],
            [AS_HELP_STRING([--with-convolution], [Choose which convolution implementation to use])],
            [use_convolution=$withval],
            [use_convolution=5])
AC_MSG_RESULT([$use_convolution])
AC_DEFINE_UNQUOTED([SPAMM_MULTIPLY_CONVOLUTE_IMPLEMENTATION], [$use_convolution], [The convolution implementation])

AC_MSG_CHECKING([whether to double-check the multiply])
AC_ARG_ENABLE([double-check],
              [AS_HELP_STRING([--enable-double-check], [Double check a few things in the multiply])],
              [use_double_check=$enableval],
              [use_double_check="no"])
AC_MSG_RESULT([$use_double_check])
if test "${use_double_check}" = "yes"; then
  AC_DEFINE([SPAMM_MULTIPLY_DOUBLE_CHECK], [1], [Double check a few things in the multiply])
fi

dnl Fix the package version in the documentation.
AC_SUBST(PACKAGE_VERSION_STRING, [${PACKAGE_VERSION}])
AC_DEFINE_UNQUOTED(SPAMM_VERSION, ["${VERSION}"], [The spammpack version.])

dnl Store the command line.
AC_DEFINE_UNQUOTED([SPAMM_CONFIGURE_ARGUMENTS], ["${ac_configure_args}"], [The command line with which spammpack was configured.])

dnl Get the git commit tag.
AC_MSG_NOTICE([sourcing ${srcdir}/git_commit])
if test -f ${srcdir}/git_commit; then
  source ${srcdir}/git_commit
else
  git_commit="unknown"
fi
AC_SUBST([GIT_COMMIT], [${git_commit}])

dnl Get the current date.
timestamp=`date +%Y-%m-%d`
AC_SUBST([TIMESTAMP], [${timestamp}])

dnl Write results.
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
Doxyfile
Makefile
interfaces/Makefile
interfaces/python/Makefile
src/Makefile
src/spamm_manual_mainpage.h
tests/Makefile
])

AC_CONFIG_COMMANDS([src/spamm_config.h], [
  echo -n "configure: creating src/spamm_config.h... "
  ./extract_conf.py config.h > src/spamm_config.h.temp
  if (diff --brief src/spamm_config.h src/spamm_config.h.temp > /dev/null 2>&1); then
    echo -n "no difference, keeping old one... "
    rm -f src/spamm_config.h.temp
  else
    echo -n "updating src/spamm_config.h... "
    mv -f src/spamm_config.h.temp src/spamm_config.h
  fi
  echo -n "removing duplicate entries from config.h... "
  sed -e "s/^#define SPAMM_/\/\/#define SPAMM_/" -i config.h
  echo "done"
])

AC_OUTPUT
