AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CLEANFILES = *~

DISTCLEANFILES = config_fortran.h

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src interfaces tests prototype

CONFIG_STATUS_DEPENDENCIES = \
  $(top_srcdir)/Doxyfile.in \
  $(top_builddir)/git_commit_tag \
  $(top_srcdir)/src/spamm_manual_mainpage.h.in

EXTRA_DIST = \
  $(DOCFILES) \
  extract_conf.py \
  update_git_commit_tag.sh \
  $(builddir)/git_commit_tag

DISTCHECK_CONFIGURE_FLAGS = --with-norm-type=double

include_HEADERS = config.h config_fortran.h

TAGS_FILES = src/* tests/*
AM_CTAGSFLAGS = --C-kinds=+l+x

#DOCFILES = latex/refman.pdf

latex/refman.pdf : docs

docs :
	doxygen
	$(MAKE) -C latex
