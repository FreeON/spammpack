SUBDIRS = libchunk src tools tests

#DISTCHECK_CONFIGURE_FLAGS = --enable-debug

EXTRA_DIST = \
	configure.LB \
	configure.mustang \
	freeon.job.template \
	freeon-submit.sh \
	moonlight.sh.template \
	moonlight-submit.sh \
	mustang.job.decay.template \
	mustang.job.diagonal.template \
	mustang.job.full.template \
	mustang.job.PP.template \
	mustang.job.water.template \
	mustang-results.sh \
	mustang-submit.sh \
	mustang-water-openmp-submit.sh \
	mustang-water-serial-submit.sh \
	openMP.job \
	spamm.bib

docs :
	$(DOXYGEN)

docs-pdf : docs
	$(MAKE) -C latex
