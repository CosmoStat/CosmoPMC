# $COSMOPMC/Makefile
# Master CosmoPMC Makefile.
# Karim Benabed, Martin Kilbinger 2010, 2011


export COSMOPMC = $(shell pwd)


include Makefile.main


ifeq ($(DOWMAP),1)
   MYCMB = cmb
else
   MYCMB =
endif


ALL = tools exec links_exec links_R links_Demo links_MC_Demo Manual
all: $(ALL)

PAR = cosmo.par cosmo_SN.par cosmo_lens.par cosmoDP.par halomodel.par cosmo_3rd.par

Makefile.host:
	@echo "Run '[python] ./configure.py' to create 'Makefile.host'"
	@exit 1

dummy:

tools: dummy Makefile.host
	cd $(TOOLS) && $(MAKE)

wrappers: dummy
	cd $(WRAPPERS) && $(MAKE)

exec: dummy tools wrappers $(MYCMB)
	cd $(EXEC) && $(MAKE)

links_exec: exec
	for target in $(EXEC_TARGETS); do \
	   ln -sf ../exec/$$target bin/$$target; \
	done

links_R:
	ln -sf ../R/plot_confidence.R bin/plot_confidence.R; \
	ln -sf ../R/sample_from_pmcsimu.R bin/sample_from_pmcsimu.R

links_Demo: 
ifeq ($(DOWMAP), 1)
	ln -sf $(COSMOPMC)/CMB/scamb Demo
	ln -sf $(MY_WMAP) Demo/WMAP
	$(COSMOPMC)/bin/config_pmc_to_max_and_fish.pl -c $(DEMO)/MC_Demo/WMAP7/config_pmc -M -f "0.043 0.25 0.087 0.96 2.4 0.7 -1.0" > $(DEMO)/config_max_WMAP7
	ln -sf $(COSMOPMC)/data/WMAP/wmap_clsz_V_l=2-10000_v4.txt Demo
endif
	for par in $(PAR); do \
	   ln -s $(COSMOPMC)/par_files/$$par Demo; \
	done

links_MC_Demo:
	# Data files
	ln -sf $(COSMOPMC)/data/BAO/bao_Reid10_A_0.35 Demo/MC_Demo/BAO/distance_A
	ln -sf $(COSMOPMC)/data/BAO/bao_BOSS12_d_z_0.57 Demo/MC_Demo/BAO/distance_d_z
	ln -sf $(COSMOPMC)/data/BAO/bao_Reid10_A_0.35 Demo/MC_Demo/COSMOS-S10+SN+BAO/
ifeq ($(DOWMAP), 1)
	ln -sf $(COSMOPMC)/data/WMAP/wmap_clsz_V_l=2-10000_v4.txt Demo/MC_Demo/WMAP$(NYR)
	ln -sf $(COSMOPMC)/CMB/scamb Demo/MC_Demo/WMAP$(NYR)
	ln -sf $(MY_WMAP) Demo/MC_Demo/WMAP$(NYR)/WMAP
endif
	for target in $(COSMOPMC)/data/lens/COSMOS/*; do \
           ln -sf $$target Demo/MC_Demo/Lensing/COSMOS-S10/`basename $$target`; \
           ln -sf $$target Demo/MC_Demo/COSMOS-S10+SN+BAO/`basename $$target`; \
	done
	for target in $(COSMOPMC)/data/lens/CFHTLenS-K13/*; do \
           ln -sf $$target Demo/MC_Demo/Lensing/CFHTLenS-K13/`basename $$target`; \
        done
	for target in $(COSMOPMC)/data/lens/CFHTLenS-FK14/*; do \
           ln -sf $$target Demo/MC_Demo/Lensing/CFHTLenS-FK14/`basename $$target`; \
        done
	for target in $(COSMOPMC)/data/HOD/berwein02_hexcl-PMC/*; do \
           ln -sf $$target Demo/MC_Demo/HOD/berwein02_hexcl-PMC/`basename $$target`; \
	done
	for target in $(COSMOPMC)/data/HOD/CFHTLS-T06/*; do \
           ln -sf $$target Demo/MC_Demo/HOD/CFHTLS-T06/`basename $$target`; \
        done
	ln -sf $(COSMOPMC)/data/Sn/Union/sne_union_marek.list Demo/MC_Demo/SN
	ln -sf $(COSMOPMC)/data/Sn/Union/sne_union_marek.list Demo/MC_Demo/COSMOS-S10+SN+BAO/
	ln -sf $(COSMOPMC)/data/WMAP_Distance_Priors/wmap7DistPrior_ML_covinv Demo/MC_Demo/WMAP_Distance_Priors

        # Parameter files
	ln -sf $(COSMOPMC)/par_files/cosmoDP.par Demo/MC_Demo/BAO/distance_A
	ln -sf $(COSMOPMC)/par_files/cosmoDP.par Demo/MC_Demo/BAO/distance_d_z
	ln -sf $(COSMOPMC)/par_files/cosmoDP.par Demo/MC_Demo/WMAP$(NYR)
	ln -sf $(COSMOPMC)/par_files/cosmo_SN.par Demo/MC_Demo/SN
	ln -sf $(COSMOPMC)/par_files/cosmo.par Demo/MC_Demo/SN
	ln -sf $(COSMOPMC)/par_files/cosmoDP.par Demo/MC_Demo/WMAP_Distance_Priors
	ln -sf $(COSMOPMC)/par_files/cosmoDP.par Demo/MC_Demo/COSMOS-S10+SN+BAO
	ln -sf $(COSMOPMC)/par_files/cosmo_SN.par Demo/MC_Demo/COSMOS-S10+SN+BAO
	ln -sf $(COSMOPMC)/par_files/cosmo.par Demo/MC_Demo/COSMOS-S10+SN+BAO
	ln -sf $(COSMOPMC)/par_files/cosmo.par Demo/MC_Demo/HOD/CFHTLS-T06
	ln -sf $(COSMOPMC)/par_files/cosmo.par Demo/MC_Demo/Lensing/CFHTLenS-K13
	ln -sf $(COSMOPMC)/par_files/cosmo.par Demo/MC_Demo/Lensing/CFHTLenS-FK14
	ln -sf $(COSMOPMC)/par_files/cosmo_lens.par Demo/MC_Demo/Lensing/CFHTLenS-K13
	ln -sf $(COSMOPMC)/par_files/cosmo_3rd.par Demo/MC_Demo/Lensing/CFHTLenS-FK14

cmb: tools
	cd $(CMB) && $(MAKE) scamb $(LIBWMAP)

manual:
	cd Manual && $(MAKE)

subdirs  = $(TOOLS) $(WRAPPERS) $(DEMO) $(EXEC)

ifeq ($(PMCLIB_EXT),0)
   subdirs += $(PMC) 
endif

ifeq ($(DOWMAP),1)
   subdirs += $(CMB)
endif


.PHONY : clean clean_links

clean:
	for dir in $(subdirs); do \
	   $(MAKE) -C $$dir clean; \
	done
	$(MAKE) clean_links
	cd Manual; $(MAKE) clean

clean_links:
	for target in $(EXEC_TARGETS); do \
	   rm -f bin/$$target; \
        done
	#rm -f sample_from_pmcsimu.R plot_confidence.R
	rm -f Demo/scamb Demo/WMAP
	rm -f Demo/MC_Demo/WMAP$(NYR)/wmap_clsz_V_l=2-10000_v4.txt Demo/MC_Demo/WMAP$(NYR)/scamb Demo/MC_Demo/WMAP$(NYR)/WMAP
	rm -f $(addprefix Demo/, $(PAR))
	for target in $(COSMOPMC)/data/lens/COSMOS/*; do \
	   rm -f Demo/MC_Demo/Lensing/COSMOS-S10/`basename $$target`; \
	   rm -f Demo/MC_Demo/COSMOS-S10+SN+BAO/`basename $$target`; \
	done
	for target in $(COSMOPMC)/data/lens/CFHTLenS-K13/*; do \
	  rm -f Demo/MC_Demo/Lensing/CFHTLenS-K13/`basename $$target`; \
	done
	for target in $(COSMOPMC)/data/lens/CFHTLenS-FK14/*; do \
           rm -f Demo/MC_Demo/Lensing/CFHTLenS-FK14/`basename $$target`; \
        done
	for target in $(COSMOPMC)/data/HOD/berwein02_hexcl-PMC/*; do \
           rm -f Demo/MC_Demo/HOD/berwein02_hexcl-PMC/`basename $$target`; \
	done
	for target in $(COSMOPMC)/data/HOD/CFHTLS-T06/*; do \
           rm -f Demo/MC_Demo/HOD/CFHTLS-T06/`basename $$target`; \
        done
	rm -f Demo/MC_Demo/SN/sne_union_marek.list Demo/MC_Demo/COSMOS-S10+SN+BAO/sne_union_marek.list
	rm -f Demo/MC_Demo/WMAP_Distance_Priors/wmap7DistPrior_ML_covinv
	rm -f Demo/MC_Demo/BAO/distance_A/cosmoDP.par Demo/MC_Demo/BAO/distance_d_z/cosmoDP.par Demo/MC_Demo/WMAP$(NYR)/cosmoDP.par
	rm -f Demo/MC_Demo/BAO/distance_d_z/bao_BOSS12_d_z_0.57 Demo/MC_Demo/BAO/distance_A/bao_Reid10_A_0.35
	rm -f Demo/MC_Demo/SN/cosmo_SN.par Demo/MC_Demo/SN/cosmo.par
	rm -f Demo/MC_Demo/WMAP_Distance_Priors/cosmoDP.par
	rm -f Demo/MC_Demo/COSMOS-S10+SN+BAO/cosmoDP.par Demo/MC_Demo/COSMOS-S10+SN+BAO/cosmo_SN.par
	rm -f Demo/MC_Demo/COSMOS-S10+SN+BAO/cosmo.par Demo/MC_Demo/COSMOS-S10+SN+BAO/bao_Reid10_A_0.35
	rm -f Demo/MC_Demo/Lensing/CFHTLenS-K13/cosmo.par
	rm -f Demo/MC_Demo/Lensing/CFHTLenS-FK14/cosmo.par
	rm -f Demo/MC_Demo/Lensing/CFHTLenS-K13/cosmo_lens.par
	rm -f Demo/MC_Demo/Lensing/CFHTLenS-FK14/cosmo_3rd.par
	rm -f Demo/MC_Demo/HOD/CFHTLS-T06/cosmo.par


