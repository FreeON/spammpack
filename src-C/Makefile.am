AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CLEANFILES = *~

DISTCLEANFILES = config_fortran.h git_commit_tag

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src interfaces tools tests

CONFIG_STATUS_DEPENDENCIES = \
  $(top_srcdir)/Doxyfile.in \
  $(top_builddir)/git_commit_tag \
  $(top_srcdir)/manual.h.in

EXTRA_DIST = \
  extract_conf.py \
  update_git_commit_tag.sh \
  $(builddir)/git_commit_tag

TAGS_FILES = src/* tests/*
AM_CTAGSFLAGS = --C-kinds=+l+x

latex/refman.pdf : docs

docs :
	doxygen
