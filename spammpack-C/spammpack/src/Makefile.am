# The version scheme used by Libtool tracks interfaces, where an interface is
# the set of exported entry points into the library. All Libtool libraries
# start with `-version-info' set to `0:0:0' -- this will be the default
# version number if you don't explicitly set it on the Libtool link command
# line. The meaning of these numbers (from left to right) is as follows:
#
# current
#     The number of the current interface exported by the library. A current
#     value of `0', means that you are calling the interface exported by this
#     library interface 0.
#
# revision
#     The implementation number of the most recent interface exported by this
#     library. In this case, a revision value of `0' means that this is the
#     first implementation of the interface.
#
#     If the next release of this library exports the same interface, but has
#     a different implementation (perhaps some bugs have been fixed), the
#     revision number will be higher, but current number will be the same. In
#     that case, when given a choice, the library with the highest revision
#     will always be used by the runtime loader.
#
# age
#     The number of previous additional interfaces supported by this library.
#     If age were `2', then this library can be linked into executables which
#     were built with a release of this library that exported the current
#     interface number, current, or any of the previous two interfaces.
#
#     By definition age must be less than or equal to current. At the outset,
#     only the first ever interface is implemented, so age can only be `0'.
#
# For later releases of a library, the `-version-info' argument needs to be
# set correctly depending on any interface changes you have made. This is
# quite straightforward when you understand what the three numbers mean:
#
#    1. If you have changed any of the sources for this library, the revision
#    number must be incremented. This is a new revision of the current
#    interface.
#
#    2. If the interface has changed, then current must be incremented, and
#    revision reset to `0'. This is the first revision of a new interface.
#
#    3. If the new interface is a superset of the previous interface (that is,
#    if the previous interface has not been broken by the changes in this new
#    release), then age must be incremented. This release is backwards
#    compatible with the previous release.
#
#    4. If the new interface has removed elements with respect to the previous
#    interface, then you have broken backward compatibility and age must be
#    reset to `0'. This release has a new, but backwards incompatible
#    interface.
#
#     For example, if the next release of the library included some new
#     commands for an existing socket protocol, you would use -version-info
#     1:0:1. This is the first revision of a new interface. This release is
#     backwards compatible with the previous release.
#
#     Later, you implement a faster way of handling part of the algorithm at
#     the core of the library, and release it with -version-info 1:1:1. This
#     is a new revision of the current interface.
#
#     Unfortunately the speed of your new implementation can only be fully
#     exploited by changing the API to access the structures at a lower level,
#     which breaks compatibility with the previous interface, so you release
#     it as -version-info 2:0:0. This release has a new, but backwards
#     incompatible interface.
#
# When deciding which numbers to change in the -version-info argument for a
# new release, you must remember that an interface change is not limited to
# the API of the library. The notion of an interface must include any method
# by which a user (code or human) can interact with the library: adding new
# builtin commands to a shell library; the format used in an output file; the
# handshake protocol required for a client connecting over a socket, and so
# on.
#
# Additionally, If you use a development model which has both a stable and an
# unstable tree being developed in parallel, for example, and you don't mind
# forcing your users to relink all of the applications which use one of your
# Libtool libraries every time you make a release, then libtool provides the
# `-release' flag to encode the project version number in the name of the
# library, See section 11.2.1 Creating Libtool Libraries with Automake. This
# can save you library compatibility problems later if you need to, say, make
# a patch release of an older revision of your library, but the library
# version number that you should use has already been taken by another earlier
# release. In this case, you could be fairly certain that library releases
# from the unstable branch will not be binary compatible with the stable
# releases, so you could make all the stable releases with `-release 1.0' and
# begin the first unstable release with `-release 1.1'.

AM_CTAGSFLAGS = --C-kinds=+l+x
TAGS_FILES = $(top_builddir)/config.h generate_SSE_assembly.py

AM_CFLAGS = $(OPENMP_CFLAGS)
AM_CPPFLAGS = $(OPENMP_CFLAGS)
AM_FCFLAGS = $(OPENMP_FCFLAGS)

ACLOCAL_AMFLAGS = -I$(top_srcdir)/m4

lib_LTLIBRARIES = libspamm.la libspammpack_fortran.la
noinst_LTLIBRARIES = libspamm_tree.la libspamm_kernel.la

if INTERNAL_LAPACK
LIBADD_INTERNAL_LAPACK = \
  $(top_builddir)/../lapack/lapack/libfreeonlapack.la \
  $(top_builddir)/../lapack/blas/libfreeonblas.la \
  $(top_builddir)/../lapack/install/libfreeonlapackextra.la
endif

libspamm_la_LDFLAGS = -version-info 0:0:0
libspamm_la_SOURCES =
libspamm_la_LIBADD  = libspamm_tree.la libspamm_kernel.la $(LIBADD_INTERNAL_LAPACK)
#libspamm_la_DEPENDENCIES = $(LIBADD_INTERNAL_LAPACK)

noinst_SCRIPTS = generate_SSE_assembly.py spammOffsets.py
noinst_PROGRAMS = print_data_sizes print_chunk_offsets

print_chunk_offsets_LDFLAGS = libspamm.la

EXTRA_DIST = spamm_list_driver.c generate_SSE_assembly.py

spammOffsets.py : print_data_sizes$(EXEEXT)
	$(builddir)/print_data_sizes$(EXEEXT) > $@

spamm_stream_kernel_SSE.S : $(srcdir)/generate_SSE_assembly.py spammOffsets.py
	PYTHONPATH=$(builddir) $(PYTHON) $(srcdir)/generate_SSE_assembly.py -N 4 --SSE 1 --name spamm_stream_kernel_SSE > $@

spamm_stream_kernel_SSE4_1.S : $(srcdir)/generate_SSE_assembly.py spammOffsets.py
	PYTHONPATH=$(builddir) $(PYTHON) $(srcdir)/generate_SSE_assembly.py -N 4 --SSE 4.1 --name spamm_stream_kernel_SSE4_1 > $@

generated_kernel_files = \
  spamm_stream_kernel_SSE.S \
  spamm_stream_kernel_SSE4_1.S

BUILT_SOURCES = $(generated_kernel_files) spammpack.mod

spammpack.mod : libspammpack_fortran.la

CLEANFILES = \
  $(generated_kernel_files) \
  *.mod

DISTCLEANFILES = \
  $(builddir)/spammOffsets.py \
  $(builddir)/spammOffsets.pyc

nodist_libspamm_kernel_la_SOURCES = $(generated_kernel_files)

libspammpack_fortran_la_SOURCES = \
  spamm_fortran_chunk.F90 \
  spamm_fortran_derived.F90 \
  spamm_fortran_interface.c \
  spamm_fortran_interface.h \
  spammpack.F90

libspamm_kernel_la_SOURCES = \
  $(libspamm_kernel_headers) \
  spamm_stream_external_sgemm.c \
  spamm_kernel.c

libspamm_kernel_headers = \
  spamm_kernel.h

libspamm_tree_la_SOURCES = \
  $(libspamm_tree_headers) \
  $(libspamm_tree_private_headers) \
  spamm_manual_mainpage.h \
  spamm_add.c \
  spamm_allocate.c \
  spamm_check.c \
  spamm_chunk.c \
  spamm_convert.c \
  spamm_copy.c \
  spamm_delete.c \
  spamm_error.c \
  spamm_get.c \
  spamm_hashtable.c \
  spamm_index.c \
  spamm_index_2D.c \
  spamm_index_3D.c \
  spamm_ipow.c \
  spamm_list.c \
  spamm_maintenance.c \
  spamm_multiply.c \
  spamm_new.c \
  spamm_print.c \
  spamm_set.c \
  spamm_sgemm.c \
  spamm_stats.c \
  spamm_timer.c \
  spamm_uint_to_bin_string.c \
  spamm_version.c

libspamm_tree_headers = \
  spamm.h \
  spamm_chunk.h \
  spamm_error.h \
  spamm_errno.h \
  spamm_general.h \
  spamm_hashtable.h \
  spamm_list.h \
  spamm_timer.h \
  spamm_types.h

libspamm_tree_private_headers = \
  spamm_types_private.h

include_HEADERS = \
  $(libspamm_kernel_headers) \
  $(libspamm_tree_headers) \
  spammpack.mod
