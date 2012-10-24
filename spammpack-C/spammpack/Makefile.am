AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CLEANFILES = *~

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src interfaces tests

CONFIG_STATUS_DEPENDENCIES = \
  $(top_srcdir)/Doxyfile.in \
  $(top_srcdir)/git_commit \
  $(top_srcdir)/src/spamm_manual_mainpage.h.in

EXTRA_DIST = latex/refman.pdf \
  extract_conf.py \
  git_commit

TAGS_FILES = src/* tests/*
AM_CTAGSFLAGS = --C-kinds=+l+x

docs :
	doxygen
	$(MAKE) -C latex

.PHONY : always
git_commit : always
	COMMIT=`git show --pretty=oneline HEAD | head -n1 | awk '{print $$1}'` && echo "git_commit=$$COMMIT" > git_commit.temp && ( diff --brief git_commit git_commit.temp > /dev/null 2>&1 ) || mv -f git_commit.temp git_commit
	@rm -f git_commit.temp
