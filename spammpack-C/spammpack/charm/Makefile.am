TAGS_FILES = config.h

DISTCHECK_CONFIGURE_FLAGS = --enable-debug

EXTRA_DIST = \
	$(INTERFACE_FILES) \
	spamm.bib \
	configure.LB \
	configure.mustang \
	freeon-submit.sh \
	freeon.job.template \
	moonlight-submit.sh \
	moonlight.sh.template \
	mustang-results.sh \
	mustang-submit.sh \
	mustang.job.decay.template \
	mustang.job.diagonal.template \
	mustang.job.full.template

AM_LDFLAGS = -language charm++ @LOADBALANCER_FLAGS@ @DEFAULT_LB_FLAGS@ @CHARMDEBUG_FLAGS@

if HAVE_PYTHON
python_PYTHON = plotPEMap.py
endif

bin_PROGRAMS = spamm BCSRInfo

INTERFACE_FILES = \
	bcsrinfo.ci \
	matrix.ci \
	messages.ci \
	multiply.ci \
	multiplyelement.ci \
	node.ci \
	spamm.ci

spamm_SOURCES = \
	bcsr.cc \
	bcsr.h \
	blas_interface.h \
	index.h \
	logger.cc \
	logger.h \
	spamm.cc \
	spamm.h \
	matrix.cc \
	matrix.h \
	messages.cc \
	messages.h \
	multiply.cc \
	multiplyelement.cc \
	multiplyelement.h \
	multiply.h \
	node.cc \
	node.h \
	timer.cc \
	timer.h \
	types.h \
	utilities.cc \
	utilities.h

$(spamm_OBJECTS) : $(BUILT_SOURCES)

BCSRInfo_SOURCES = bcsrinfo.cc bcsr.cc bcsr.h logger.h logger.cc

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

BUILT_SOURCES = $(foreach i, $(INTERFACE_FILES), $(basename $i).decl.h $(basename $i).def.h)

CLEANFILES = $(BUILT_SOURCES) charmrun
