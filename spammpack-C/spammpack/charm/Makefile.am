if HAVE_CHARM

include $(top_srcdir)/Makefile.charm

bin_PROGRAMS = multiply

INTERFACE_FILES = \
	matrix.ci \
	matrixnode.ci \
	messages.ci \
	multiply.ci

HEADER_FILES = \
	algebra.h \
	matrix.h \
	matrixnode.h \
	messages.h \
	timer.h \
	utilities.h

multiply_SOURCES = \
	$(INTERFACE_FILES) \
	$(HEADER_FILES) \
	algebra.cc \
	matrix.cc \
	matrixnode.cc \
	messages.cc \
	multiply.cc \
	timer.cc \
	utilities.cc

endif
