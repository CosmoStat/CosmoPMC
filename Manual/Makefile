# Makefile for manual
# COSMOPMC/manual/Makefile
# Martin Kilbinger 2010

include $(COSMOPMC)/Makefile.main


source = manual
DD     = $(source)-`date +"%Y%m%d"`.tgz
verb   = \
	 add_comment_to_ps.pl \
	 add_ded.pl \
	 add_deduced_cosebis \
	 add_deduced_halomodel \
         add_par_from_prior.pl \
	 add_par_to_mvdens.pl \
	 add_pmc_proposal \
	 all_vs_all.pl \
	 allps2tex.pl \
	 ascii2sample.pl \
         bayes_factor.pl \
         cl_one_sided \
         config_pmc_to_max_and_fish.pl \
         corr_coeff.sh \
         cosmo_mcmc \
         cosmo_pmc \
         cosmo_pmc.pl \
         diag_mvdens.pl \
         essential_cosmo_pmc_run.pl \
         evidence.pl \
         evidence_list.pl \
         fisher_to_meanvar.pl \
         get_spar.pl \
         go_fishing \
	 haloplot \
	 histogram.pl \
         histograms_sample \
         importance_sample \
         max_post \
         mean2eps.pl \
         meanvar2tab.pl \
         meanvar_sample \
         neff_proposal.pl \
         newdir_pmc.sh \
	 plot_confidence.R \
         plot_contour2d.pl \
	 pmcsim_cat.pl \
         proposal_mean.pl \
         proposal_var.pl \
         remap.sh \
         sample2fixpar.pl \
         sample_from_pmcsimu.R \
         tab2tex.pl \
         test_suite_cosmo_pmc.pl \
	 xi2xi.pl \


CosmoPMCVersion = `perl -F\" -ane 'print $$F[1] if /COSMO_PMC_V/' $(COSMOPMC)/tools/include/stdnames.h`


ifeq ($(STATUS), "RELEASED")
   finish = manual
else
   finish = bib
endif


all: version $(verb) $(finish) shortguide clear_version

allverb: version $(verb) clear_version


version:
	@ echo $(CosmoPMCVersion) > CosmoPMCVersion.tmp

clear_version:
	@ rm -rf CosmoPMCVersion.tmp

$(source):
	$(MAKE) version
	pdflatex -enable-write18 $(source).tex
	#pdflatex -enable-write18 $(source).tex
	#pdflatex -enable-write18 $(source).tex
	$(MAKE) clear_version

bib:
	$(MAKE) version
	$(MAKE) $(source)
	bibtex $(source)
	$(MAKE) $(source)
	$(MAKE) $(source)
	$(MAKE) clear_version

quick:
	pdflatex -enable-write18 $(source).tex

shortguide:
	$(MAKE) version
	pdflatex shortguide.tex
	pdflatex shortguide.tex
	pdflatex shortguide.tex
	$(MAKE) clear_version	

$(verb):
	./verbatim-usage.pl ../bin/$@ > usage_$@.tex 2>&1

PHONY: clean clean_pdf backup arxiv

clean:
	rm -f *~ $(source).out $(source).blg $(source).bbl *.aux *.log $(source).dvi $(source).toc $(source).out \
	  $(source).brf $(source).ps $(source).toc usage_*.tex shortguide.out
	$(MAKE) clear_version

clean_pdf: clean
	rm -f manual.pdf

arxiv:
	tar czf upload.tgz manual.tex manual.bbl usage*.tex texititel.sty preliminaries.tex abbr-journals.tex \
	 download_install.tex quick_ref.tex CosmoPMCVersion.tmp chain_flow.pdf cosmo_pmc_flow.pdf

backup:
	tar czf $(DD) --exclude *.tgz *
#	h: dereference links
	mutt -s "backup $(source)" -a $(DD) $(MAILID) < /dev/null
	cp $(DD) $(HOME)/sicher
