TAGS_FILES = config.h

EXTRA_DIST = $(spamm_INTERFACEFILES)

additional_linkerflags = @LOADBALANCER_FLAGS@ @DEFAULT_LB_FLAGS@ @CHARMDEBUG_FLAGS@

check_SCRIPTS = \
	$(top_builddir)/spamm_test_add.sh \
	$(top_builddir)/spamm_test_multiply.sh \
	$(top_builddir)/spamm_test_trace.sh
TESTS = $(check_SCRIPTS)

bin_PROGRAMS = spamm

spamm_LDFLAGS = -language charm++ $(additional_linkerflags)

spamm_INTERFACEFILES = \
	matrix.ci \
	messages.ci \
	multiply.ci \
	multiplyelement.ci \
	node.ci \
	spamm.ci

spamm_SOURCES = \
	bcsr.cc \
	bcsr.h \
	index.h \
	lapack_interface.h \
	logger.cc \
	logger.h \
	matrix.cc \
	matrix.h \
	messages.cc \
	messages.h \
	multiply.cc \
	multiply.h \
	multiplyelement.cc \
	multiplyelement.h \
	node.cc \
	node.h \
	spamm.cc \
	spamm.h \
	timer.cc \
	timer.h \
	types.h \
	utilities.cc \
	utilities.h

spamm_CHARMHEADERS = $(foreach i, $(spamm_INTERFACEFILES), $(basename $i).decl.h $(basename $i).def.h)

$(spamm_OBJECTS) : $(spamm_CHARMHEADERS)

# Matrices are blocked into 32z32 at the leaf level.
N = 4096
B = 128

RSH = rsh

run : spamm$(EXEEXT)
	for i in 1 2 4; do \
		./charmrun +p$$i ++remote-shell $(RSH) spamm -N $(N) -b $(B); \
	done

run-local : spamm$(EXEEXT)
	for i in 1 2 4; do \
		./spamm +p$$i -N $(N) -b $(B); \
	done

run-LB : spamm$(EXEEXT)
	./charmrun +p4 spamm -N $(N) -b $(B) --verify +LBDebug

valgrind : spamm$(EXEEXT)
	valgrind --log-file=valgrind.%p.log \
		--track-origins=yes \
		--leak-check=full \
		--show-reachable=yes \
		./$< -N 1

docs :
	$(DOXYGEN)

docs-pdf : docs
	$(MAKE) -C latex

%.decl.h %.def.h : %.ci
	$(AM_V_CXX)$(CHARMC) $(DEFAULT_INCLUDES) -preprocess $<

CLEANFILES = \
	$(spamm_CHARMHEADERS) \
	charmrun

DISTCLEANFILES = $(TESTS)
