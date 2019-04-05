
# Makefile for COSMOPMC executables
# COSMOPMC/exec/Makefile
# Karim Benabed, Martin Kilbinger 2007, 2008

include $(COSMOPMC)/Makefile.main


# MPI support for cflags
IFLAGS_MPI = $(INC_MPI)
# Non-standard MPI libraries (LIBS_MPI) should be set in $COSMOPMC/Makefiles/Makefile.host
LFLAGS_MPI = $(LDIRS_MPI) $(LIBS_MPI)

OBJS_PAR = parabox.o
OBJS_MC  = montecarlo2.o nhist.o mcmc.o tools.o pmc.o distribution.o
CFLAGS  += -I$(TOOLS)/include $(INC_PMC)

# Libraries for WMAP
ifeq ($(DOWMAP),1)
   FLIBS    = -L$(FLIBPATH) -Wl,-rpath,$(IFCOREDIR) -L$(IFCOREDIR) $(LIBF90)
   WMAPLIBS = -L. -L$(MY_WMAP) -L$(CMB) -lwmap$(NYR)_wrap -lwmap$(NYR) -L$(CFITSIODIR) -lcfitsio
   LAPACK   = -Wl,-rpath,$(LAPACKDIR) -L$(LAPACKDIR) $(LIBLAPACK) \
              -Wl,-rpath,$(GUIDEDIR) -L$(GUIDEDIR) # -lmkl_em64t
   thescamb =
else
   FLIBS    =
   WMAPLIBS =
   LAPACK   =
   LIBWMAP  =
   thescamb =
endif


CWD          = $(shell pwd)
MKCFLAGS     = $(INC_PMC) -I$(WRAPPERS)/include -I$(TOOLS)/include -I$(CWD)  $(IINCDIRS) \
               -I$(FFTW)/include -I$(GSL)/include $(WARNINGS) $(ARCHFLAGS) $(IFLAGS_MPI) $(FPIC) $(DEBUG) $(OPTIM)
MKLFLAGS    += -L$(WRAPPERS) -L$(TOOLS) $(LLIBDIRS) -Wl,-rpath,$(FFTW)/lib -L$(FFTW)/lib -Wl,-rpath,$(GSL)/lib \
               -L$(GSL)/lib -lparam -lwr $(LLIBS) \
	       -lerrorio -lpar -lconf -lmvdens \
               -lgsl -lgslcblas -lfftw3 -lm $(SHARED) $(LAPACK) $(FLIBS) $(WMAPLIBS) $(FPIC) $(DEBUG) \
	       $(OPTIM)

# MKDEBUG: On sap clusters astro?, dl has to come at the end
OPTIM += -ldl

ifeq ($(DOWMAP),1)
   MKCFLAGS    += -DDOWMAP
endif

PMCLIBFLAG    = $(LIB_PMC) -lpmc -ldl
PMCMPILIBFLAG = $(PMCLIBFLAG) -lpmc_mpi

# Add topolike linking flags
OPTIM := $(OPTIM) $(MKLFLAGS_LAST)


# EXEC_TARGETS are defined in Makefile.main

all:	$(EXEC_TARGETS)


cosmo_mcmc: cosmo_mcmc.o mkmax.o exec_helper.o $(thescamb)
	$(MPICC) -o $@ $@.o mkmax.o exec_helper.o $(MKLFLAGS) $(PMCLIBFLAG) $(OPTIM)

cosmo_pmc: cosmo_pmc.o exec_helper.o $(thescamb)
	$(MPICC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(LFLAGS_MPI) $(PMCMPILIBFLAG) $(OPTIM)

max_post: max_post.o mkmax.o exec_helper.o $(thescamb)
	$(CC) -o $@ $@.o mkmax.o exec_helper.o $(MKLFLAGS) $(PMCLIBFLAG) $(OPTIM)

go_fishing: go_fishing.o mkmax.o $(thescamb)
	$(MPICC) -o $@ $@.o mkmax.o $(MKLFLAGS) $(LFLAGS_MPI) $(PMCMPILIBFLAG) $(OPTIM)

xcorr_post: xcorr_post.o exec_helper.o $(thescamb)
	$(CC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(PMCLIBFLAG)

wmap3_post: wmap3_post.o exec_helper.o $(thescamb)
	$(CC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(PMCLIBFLAG)

wmap5_post: wmap5_post.o exec_helper.o $(MKLIS) $(thescamb)
	$(CC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(PMCLIBFLAG)

meanvar_sample: meanvar_sample.o exec_helper.o $(LIBERROR) $(LIBOMATH)
	$(MPICC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(LFLAGS_MPI) $(PMCLIBFLAG)

importance_sample: importance_sample.o exec_helper.o $(LIBERROR) $(LIBOMATH)
	$(MPICC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(LFLAGS_MPI) $(PMCMPILIBFLAG)

histograms_sample: histograms_sample.o exec_helper.o 
	$(MPICC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(LFLAGS_MPI) $(PMCMPILIBFLAG)

sample_cmb: sample_cmb.o exec_helper.o 
	$(MPICC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(LFLAGS_MPI) $(PMCMPILIBFLAG)

meanvar_mixmvdens: meanvar_mixmvdens.o exec_helper.o
	$(MPICC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(LFLAGS_MPI) $(PMCMPILIBFLAG)

add_pmc_proposal: add_pmc_proposal.o exec_helper.o
	$(MPICC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(LFLAGS_MPI) $(PMCMPILIBFLAG)

cl_one_sided : cl_one_sided.o exec_helper.o
	$(MPICC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(LFLAGS_MPI) $(PMCMPILIBFLAG)

add_deduced_halomodel: add_deduced_halomodel.o exec_helper.o
	$(MPICC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(LFLAGS_MPI) $(PMCMPILIBFLAG) $(PMCLIBFLAG)

add_deduced_cosebis: add_deduced_cosebis.o exec_helper.o
	$(MPICC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(LFLAGS_MPI) $(PMCMPILIBFLAG)

getHODModel: getHODModel.o  exec_helper.o
	$(CC) -o getHODModel getHODModel.o  exec_helper.o $(MKLFLAGS) $(PMCLIBFLAG) $(OPTIM)

sample_GalCorr: sample_GalCorr.o exec_helper.o
	$(MPICC) -o $@ $@.o exec_helper.o $(MKLFLAGS) $(LFLAGS_MPI) $(PMCMPILIBFLAG)

%.o:    %.c Makefile
	$(CC) -c -o $@ $< $(MKCFLAGS)

$(LIBERROR) $(LIBMVDENS):
ifeq ($(PMCLIB_EXT),0)
	cd $(TOOLS) && $(MAKE) $@
endif

$(LIBCONF) $(LIBPAR):
	cd $(TOOLS) && $(MAKE) $@

$(LIBWRAP) $(LIBPARAM):
	cd $(WRAPPERS) && $(MAKE) $@

# Add rule for new module here:
# $(LIBNEWMODULE):
#	cd $(NEWMODULE) && $(MAKE) $@	

# New: libwmap7_wrap.a is for some reason not properly dealt with the rule $(LIBWMAP)
thescamb:
	cd $(CMB) && $(MAKE) $(LIBWMAP) scamb

%.h:


PHONY: clean

clean:
	rm -f $(EXEC_TARGETS) *.o


