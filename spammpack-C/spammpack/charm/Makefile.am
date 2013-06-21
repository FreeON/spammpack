if HAVE_CHARM

include $(top_srcdir)/Makefile.charm

bin_PROGRAMS = multiply

INTERFACE_FILES = \
	matrix.ci \
	matrixnode.ci \
	messages.ci \
	multiply.ci

HEADER_FILES = \
	matrix.h \
	matrixnode.h \
	messages.h \
	utilities.h

multiply_SOURCES = \
	$(INTERFACE_FILES) \
	$(HEADER_FILES) \
	matrix.cc \
	matrixnode.cc \
	messages.cc \
	multiply.cc \
	utilities.cc

endif
