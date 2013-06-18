if HAVE_CHARM

include $(top_srcdir)/Makefile.charm

bin_PROGRAMS = spammmultiply

INTERFACE_FILES = \
	spammmatrix.ci \
	spammmultiply.ci

HEADER_FILES = \
	spammmatrix.h

spammmultiply_SOURCES = \
	$(INTERFACE_FILES) \
	$(HEADER_FILES) \
	spammmatrix.cc \
	spammmultiply.cc

endif
