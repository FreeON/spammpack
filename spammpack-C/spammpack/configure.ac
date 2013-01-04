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

dnl Get pointer int type.
AC_TYPE_INTPTR_T

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
else
  AC_MSG_NOTICE([will use external sgemm_()])
  AC_DEFINE([ADD_SGEMM_EXTERNAL_DECLARATION], [1], [Add declaration for sgemm_() external.])
  AC_DEFINE([SGEMM], [sgemm_], [External sgemm_().])
fi

dnl Check for dgemm_().
AC_SEARCH_LIBS([dgemm_], [],
               [have_dgemm="yes"],
               [have_dgemm="no"])

if test "${have_dgemm}" = "no"; then
  AC_MSG_NOTICE([no external dgemm_(), will use internal dgemm()])
  AC_DEFINE([DGEMM], [spamm_dgemm], [SpAMM internal dgemm.])
else
  AC_MSG_NOTICE([will use external dgemm_()])
  AC_DEFINE([ADD_DGEMM_EXTERNAL_DECLARATION], [1], [Add declaration for dgemm_() external.])
  AC_DEFINE([DGEMM], [dgemm_], [External dgemm_().])
fi

dnl Check for sgeev_().
AC_SEARCH_LIBS([sgeev_], [],
               [have_sgeev="yes"],
               [have_sgeev="no"])

if test "${have_sgeev}" = "no"; then
  AC_MSG_NOTICE([no external sgeev_(), will not use sgeev()])
  AC_DEFINE([SGEEV], [//sgeev_], [not using sgeev_().])
else
  AC_MSG_NOTICE([will use external sgeev_()])
  AC_DEFINE([ADD_DGEMM_EXTERNAL_DECLARATION], [1], [Add declaration for sgeev_() external.])
  AC_DEFINE([SGEEV], [sgeev_], [External sgeev_().])
fi

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

#AC_DEFINE([SPAMM_MULTIPLY_DEBUG], [1], [Extra debugging output in the multiply.])
AC_DEFINE([RUN_ASSEMBLY_KERNEL], [1], [Run the assembly kernel.])

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
AC_DEFINE_UNQUOTED([SPAMM_TIMESTAMP], ["${timestamp}"], [The build date.])

dnl Write headers.
AC_CONFIG_HEADERS([config.h])

dnl Write Fortran friendly version of header file.
AC_CONFIG_COMMANDS([config_fortran.h], [
  echo -n "configure: creating config_fortran.h... "
  ${srcdir}/extract_conf.py config.h > config_fortran.h.temp
  if (diff --brief config_fortran.h config_fortran.h.temp > /dev/null 2>&1); then
    echo "no difference, keeping old one"
    rm -f config_fortran.h.temp
  else
    echo "updating config_fortran.h"
    mv -f config_fortran.h.temp config_fortran.h
  fi
])

dnl Write makefiles.
AC_CONFIG_FILES([
Doxyfile
Makefile
interfaces/Makefile
interfaces/python/Makefile
src/Makefile
src/spamm_manual_mainpage.h
tests/Makefile
])

AC_OUTPUT
