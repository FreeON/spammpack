if HAVE_CHARM

include $(top_srcdir)/Makefile.charm

bin_PROGRAMS = spammmultiply

INTERFACE_FILES = \
	matrixnode.ci \
	messages.ci \
	spammmatrix.ci \
	spammmultiply.ci

HEADER_FILES = \
	matrixnode.h \
	messages.h \
	spammmatrix.h

spammmultiply_SOURCES = \
	$(INTERFACE_FILES) \
	$(HEADER_FILES) \
	matrixnode.cc \
	messages.cc \
	spammmatrix.cc \
	spammmultiply.cc

endif
