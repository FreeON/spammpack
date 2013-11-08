TAGS_FILES = $(top_builddir)/config.h

EXTRA_DIST = $(spamm_charm_INTERFACEFILES)

additional_linkerflags = @LOADBALANCER_FLAGS@ @DEFAULT_LB_FLAGS@ @CHARMDEBUG_FLAGS@

check_SCRIPTS = \
	$(top_builddir)/spamm_test_add.sh \
	$(top_builddir)/spamm_test_multiply.sh \
	$(top_builddir)/spamm_test_trace.sh
TESTS = $(check_SCRIPTS)

bin_PROGRAMS = spamm_charm

spamm_charm_LDFLAGS = -language charm++ $(additional_linkerflags)

spamm_charm_INTERFACEFILES = \
	matrix.ci \
	messages.ci \
	multiply.ci \
	multiplyelement.ci \
	node.ci \
	spamm_charm.ci

spamm_charm_SOURCES = \
	bcsr.cc \
	bcsr.h \
	block.cc \
	block.h \
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
	spamm_charm.cc \
	spamm_charm.h \
	timer.cc \
	timer.h \
	types.h \
	utilities.cc \
	utilities.h

spamm_charm_CHARMHEADERS = $(foreach i, $(spamm_charm_INTERFACEFILES), $(basename $i).decl.h $(basename $i).def.h)

$(spamm_charm_OBJECTS) : $(spamm_charm_CHARMHEADERS)

# Matrices are blocked into 32z32 at the leaf level.
N = 4096
B = 128

RSH = rsh

run : spamm_charm$(EXEEXT)
	for i in 1 2 4; do \
		./charmrun +p$$i ++remote-shell $(RSH) spamm_charm -N $(N) -b $(B); \
	done

run-local : spamm_charm$(EXEEXT)
	for i in 1 2 4; do \
		./spamm_charm +p$$i -N $(N) -b $(B); \
	done

run-LB : spamm_charm$(EXEEXT)
	./charmrun +p4 spamm_charm -N $(N) -b $(B) --verify +LBDebug

valgrind : spamm_charm$(EXEEXT)
	valgrind --log-file=valgrind.%p.log \
		--track-origins=yes \
		--leak-check=full \
		--show-reachable=yes \
		./$< -N 1

%.decl.h %.def.h : %.ci
	$(AM_V_CXX)$(CHARMC) $(DEFAULT_INCLUDES) -preprocess $<

CLEANFILES = \
	$(spamm_charm_CHARMHEADERS) \
	charmrun

DISTCLEANFILES = $(TESTS)
