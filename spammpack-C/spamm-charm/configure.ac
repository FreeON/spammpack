AC_INIT([spamm-charm],
        [1],
        [nicolas.bock@freeon.org],
        [],
        [http://www.freeon.org/])
AM_INIT_AUTOMAKE([dist-bzip2 no-dist-gzip silent-rules])
AM_SILENT_RULES([yes])

# Find doxygen and related programs.
AC_ARG_VAR([DOXYGEN], [The doxygen executable (for building the documentation)])
AC_CHECK_PROG([DOXYGEN], [doxygen], [doxygen])

AC_CHECK_PROG([DOT], [dot], [dot])
if test -n ${DOT}; then
  have_dot="yes"
else
  have_dot="no"
fi
AC_SUBST([PROJECT_NAME], [${PACKAGE_NAME}])
AC_SUBST([PROJECT_NUMBER], [${PACKAGE_VERSION}])
AC_SUBST([HAVE_DOT], [${have_dot}])

# The doxygen configuration file is generated.
AC_SUBST([CONFIG_STATUS_DEPENDENCIES], ['$(top_srcdir)/Doxyfile.in'])

# Find a C compiler.
AC_PROG_CC_STDC

# Find OpenMP for C.
AC_LANG(C)
AC_OPENMP

# Initialize libtool.
LT_INIT

# Find a python interpreter.
AM_PATH_PYTHON([2.7], [], [:])
AM_CONDITIONAL([HAVE_PYTHON], [test "${PYTHON}" != ":"])

# Find the real C++ compiler.
AC_PROG_CXX
AC_LANG(C++)
AC_OPENMP

# Test for spammpack.
#AC_CHECK_HEADERS([spamm.h], [], [AC_MSG_FAILURE([Cannot find spammpack header files])])
#AC_MSG_CHECKING([whether I can link against spammpack])
#AC_LINK_IFELSE([AC_LANG_PROGRAM([#include <spamm.h>],
#  [  unsigned int N[[]] = { 16, 16 };
#  unsigned int N_lower[[]] = { 0, 0 };
#  unsigned int N_upper[[]] = { 16, 16 };
#  spamm_chunk_t *M = spamm_new_chunk(2, 1, N, N_lower, N_upper);])],
#  [AC_MSG_RESULT([yes])],
#  [AC_MSG_FAILURE([Cannot link against spammpack library])])

AC_SEARCH_LIBS([clock_gettime], [rt],
  [],
  [AC_MSG_FAILURE([need clock_gettime()])])

# Check for supported clocks. Thanks to Jonathan Lifflander
# <jliffl2@illinois.edu> for suggesting this test.
for CLOCKTYPE in \
  CLOCK_MONOTONIC_RAW \
  CLOCK_MONOTONIC \
  CLOCK_REALTIME; do
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#include <time.h>]],
    [[struct timespec p; clock_gettime(${CLOCKTYPE}, &p);]])],
    [break], [CLOCKTYPE="unset"])
done
if test "${CLOCKTYPE}" = "unset"; then
  AC_MSG_FAILURE([no known clock type])
fi
AC_DEFINE_UNQUOTED([CLOCKTYPE], [${CLOCKTYPE}], [The clock type to use with clock_gettime()])

AC_ARG_ENABLE([lapack],
  [AS_HELP_STRING([--enable-lapack],
  [Enable external linkage to blas/lapack])],
  [use_external_lapack=$enableval],
  [use_external_lapack="no"])

if test "${use_external_lapack}" = "yes"; then
  # Find dgemm() for matrix block products.
  for DGEMM in dgemm dgemm_ dgemm__; do
    AC_SEARCH_LIBS([${DGEMM}], [blas], [break], [DGEMM="unset"])
  done
  if test "${DGEMM}" != "unset"; then
    AC_DEFINE_UNQUOTED([DGEMM], [${DGEMM}], [The call to BLAS dgemm() including
                        all necessary underscores.])
  fi

  # Find dsyev() for spectral bounds.
  for DSYEV in dsyev dsyev_ dsyev__; do
    AC_SEARCH_LIBS([${DSYEV}], [lapack], [break], [DSYEV="unset"])
  done
  if test "${DSYEV}" != "unset"; then
    AC_DEFINE_UNQUOTED([DSYEV], [${DSYEV}], [The call to dsyev() including all
                        necessary underscores.])
  fi
fi

# Find the Charm++ compiler wrapper. We first reset the cached values autoconf
# found for the C++ compiler and then test charmc. If we don't, autoconf will
# not actually do anything, since it thinks that we have already checked.
unset ac_cv_prog_ac_ct_CXX
unset ac_cv_cxx_compiler_gnu
unset ac_cv_prog_cxx_g
unset am_cv_CXX_dependencies_compiler_type
AC_ARG_VAR([CHARMC], [The Charm++ compiler wrapper])
if test -z "${CHARMC}"; then
  CHARMC="charmc"
fi
CXX=${CHARMC}
AC_PROG_CXX

AC_LANG(C++)

# Test whether the charm compiler works.
AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <charm++.h>],
  [CkExit();])],
  [],
  [AC_MSG_FAILURE([can not compile with the charm compiler: ${CHARMC}])])

# Enable load balancers.
AC_ARG_ENABLE([load-balance],
  [AS_HELP_STRING([--enable-load-balance],
  [enable a common set of load balancers during compilation])],
  [enable_LB=$enableval],
  [enable_LB="no"])
if test "${enable_LB}" = "yes"; then
  AC_SUBST([LOADBALANCER_FLAGS], ["-module CommonLBs"])
fi

# Enable the default load balancer.
AC_ARG_ENABLE([default-LB],
  [AS_HELP_STRING([--enable-default-LB],
  [set the default load balancer (default = GreedyCommLB)])],
  [default_LB=$enableval], [])
if test -n "${default_LB}"; then
  if test "${default_LB}" = "yes"; then
    default_LB="GreedyCommLB"
  fi
  AC_SUBST([DEFAULT_LB_FLAGS], ["-balancer ${default_LB}"])
fi

# Enable charmdebug.
AC_ARG_ENABLE([charmdebug],
  [AS_HELP_STRING([--enable-charmdebug],
  [enable the use of charmdebug on the executable, use with CXXFLAGS="-O0 -g"])],
  [enable_charmdebug=$enableval],
  [enable_charmdebug="no"])
if test "${enable_charmdebug}" = "yes"; then
  AC_SUBST([CHARMDEBUG_FLAGS], ["-memory charmdebug"])
fi

# Enable debugging.
AC_MSG_CHECKING([whether to turn on full debugging output])
AC_ARG_ENABLE([debug],
  [AS_HELP_STRING([--enable-debug],
  [enable a lot of extra output])],
  [enable_debug=$enableval],
  [enable_debug="no"])
if test "${enable_debug}" = "yes"; then
  AC_DEFINE([DEBUG_OUTPUT], [1], [Enable a lot of extra output])
  AC_MSG_RESULT([yes])
else
  AC_MSG_RESULT([no])
fi

# Print some matrices for debugging.
AC_ARG_ENABLE([print-matrices],
  [AS_HELP_STRING([--enable-print-matrices],
  [print matrices for debugging])],
  [enable_print_matrices=$enableval],
  [enable_print_matrices="no"])
if test "${enable_print_matrices}" = "yes"; then
  AC_DEFINE([PRINT_MATRICES], [1], [Print matrices for debugging])
fi

# Enable pruning of the convolution. This is done by calling ckDestroy() on
# convolution elements.
AC_ARG_ENABLE([pruning],
  [AS_HELP_STRING([--enable-pruning],
  [enable pruning of convolution space])],
  [enable_pruning=$enableval],
  [enable_pruning="no"])
if test "${enable_pruning}" = "yes"; then
  AC_DEFINE([PRUNE_CONVOLUTION], [1], [Enable pruning of convolution space])
fi

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([
Makefile
src/Makefile
tools/Makefile
tests/Makefile
Doxyfile])

AC_CONFIG_COMMANDS([tests/test_sp2.sh],
[cat > tests/test_sp2.sh << END_OF_SCRIPT
#!/bin/bash
\${top_srcdir}/tools/sp2.py --Ne 2 random.mm || exit 1
END_OF_SCRIPT
chmod 755 tests/test_sp2.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_SP2_full.sh],
[cat > tests/spamm_test_SP2_full.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation SP2 -N 7 --Ne 2 --block 5 --type full --iterations 40 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_SP2_full.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_add_full.sh],
[cat > tests/spamm_test_add_full.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation add -N 21 --block 5 --type full --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_add_full.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_add_decay.sh],
[cat > tests/spamm_test_add_decay.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation add -N 21 --block 5 --type decay --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_add_decay.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_add_diagonal.sh],
[cat > tests/spamm_test_add_diagonal.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation add -N 21 --block 5 --type diagonal --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_add_diagonal.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_multiply_full.sh],
[cat > tests/spamm_test_multiply_full.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation multiply -N 21 --block 5 --type full --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_multiply_full.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_multiply_decay.sh],
[cat > tests/spamm_test_multiply_decay.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation multiply -N 21 --block 5 --type decay --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_multiply_decay.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_multiply_diagonal.sh],
[cat > tests/spamm_test_multiply_diagonal.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation multiply -N 21 --block 5 --type diagonal --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_multiply_diagonal.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_scale_full.sh],
[cat > tests/spamm_test_scale_full.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation scale -N 21 --block 5 --type full --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_scale_full.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_scale_decay.sh],
[cat > tests/spamm_test_scale_decay.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation scale -N 21 --block 5 --type decay --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_scale_decay.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_scale_diagonal.sh],
[cat > tests/spamm_test_scale_diagonal.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation scale -N 21 --block 5 --type diagonal --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_scale_diagonal.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_trace_full.sh],
[cat > tests/spamm_test_trace_full.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation trace -N 21 --block 5 --type full --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_trace_full.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_trace_decay.sh],
[cat > tests/spamm_test_trace_decay.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation trace -N 21 --block 5 --type decay --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_trace_decay.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_trace_diagonal.sh],
[cat > tests/spamm_test_trace_diagonal.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation trace -N 21 --block 5 --type diagonal --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_trace_diagonal.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_addIdentity_full.sh],
[cat > tests/spamm_test_addIdentity_full.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation addIdentity -N 21 --block 5 --type full --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_addIdentity_full.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_addIdentity_decay.sh],
[cat > tests/spamm_test_addIdentity_decay.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation addIdentity -N 21 --block 5 --type decay --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_addIdentity_decay.sh
])

AC_CONFIG_COMMANDS([tests/spamm_test_addIdentity_diagonal.sh],
[cat > tests/spamm_test_addIdentity_diagonal.sh << END_OF_SCRIPT
#!/bin/bash
\${top_builddir}/src/charmrun +p3 \${top_builddir}/src/spamm-charm --operation addIdentity -N 21 --block 5 --type diagonal --iterations 3 --load-balance --verify || exit 1
END_OF_SCRIPT
chmod 755 tests/spamm_test_addIdentity_diagonal.sh
])

AC_CONFIG_COMMANDS([tests/random.mm],
[cat > tests/random.mm << END_OF_SCRIPT
%%MatrixMarket matrix coordinate double general
% Random matrix
3 3 9
1 1 0.18501521481190797
1 2 1.4468044758617939
1 3 1.1553404778926066
2 1 1.4468044758617939
2 2 1.9687648266446316
2 3 1.2168112202078234
3 1 1.1553404778926066
3 2 1.2168112202078234
3 3 0.97822861428160168
END_OF_SCRIPT
])

AC_OUTPUT
