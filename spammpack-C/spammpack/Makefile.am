AUTOMAKE_OPTIONS = dist-bzip2 no-dist-gzip

CLEANFILES = *~

DISTCLEANFILES = config_fortran.h

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src interfaces tests

CONFIG_STATUS_DEPENDENCIES = \
  $(top_srcdir)/Doxyfile.in \
  $(top_srcdir)/git_commit \
  $(top_srcdir)/src/spamm_manual_mainpage.h.in

EXTRA_DIST = \
  $(DOCFILES) \
  extract_conf.py \
  git_commit

include_HEADERS = config.h config_fortran.h

TAGS_FILES = src/* tests/*
AM_CTAGSFLAGS = --C-kinds=+l+x

#DOCFILES = latex/refman.pdf

latex/refman.pdf : docs

docs :
	doxygen
	$(MAKE) -C latex

.PHONY : always
git_commit : always
	@git --version > /dev/null 2>&1; if test $$? -eq 0; then \
		echo "updating git_commit"; \
		COMMIT=`git show --pretty=oneline HEAD | head -n1 | awk '{print $$1}'` && \
		echo "git_commit=$$COMMIT" > git_commit.temp && \
		( diff --brief git_commit git_commit.temp > /dev/null 2>&1 ) || mv -f git_commit.temp git_commit; \
	else \
		echo "not updating git_commit"; \
	fi
	@rm -f git_commit.temp
