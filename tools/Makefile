# Makefile for tools
# COSMOPMC/tools/Makefile
# Martin Kilbinger 2008

include $(COSMOPMC)/Makefile.main

cflags        = -Iinclude $(INC_PMC) -I$(FFTW)/include -I$(GSL)/include $(DEBUG) $(OPTIM) $(WARNINGS) $(ARCHFLAGS) $(FPIC)
lflags        = -L$(FFTW)/lib -L$(GSL)/lib -lm -lgsl -lgslcblas -lfftw3 $(FPIC) $(OPTIM)

libconfObj    = config.o nhist.o
libparObj     = par.o

ifeq ($(PMCLIB_EXT),0)
   libmvdensObj  = maths.o mvdens.o
   liberrorObj   = errorlist.o io.o
else
   libmvdenshObj  =
   liberrorObj   =
endif

vpath %.h ../tools/include
vpath %.c ../tools/src


ALL = libconf.$(EXT) libpar.$(EXT)
ifeq ($(PMCLIB_EXT),0)
   ALL += liberrorio.$(EXT) libmvdens.$(EXT)
endif

all : check_variables $(ALL)

ifeq ("$(PMC)", "")
check_variables:
	@echo "Error: Environment variable 'PMCLIB' is not set"
	@exit 1
else
check_variables:
endif


liberrorio.$(EXT_DYN): $(liberrorObj)
	$(LD) $(SHARED) -o $ $@ $^

liberrorio.a : $(liberrorObj)
	$(AR) -r $@ $^

libmvdens.$(EXT_DYN): $(libmvdensObj)
	$(LD) $(SHARED) -o $ $@ $^ $(lflags)

libmvdens.a: $(libmvdensObj)
	$(AR) -r $@ $^

libpar.$(EXT_DYN): $(libparObj)
	$(LD) $(SHARED) -o $ $@ $^ $(lflags)

libpar.a : $(libparObj)
	$(AR) -r $@ $^

libconf.$(EXT_DYN): $(libconfObj)
	$(LD) $(SHARED) -o $ $@ $^ $(lflags)

libconf.a : $(libconfObj)
	$(AR) -r $@ $^

%.o:	%.c %.h Makefile
	$(CC) $< -c $(cflags) -o $@

.h:


PHONY: clean tar

clean:
	rm -rf *~ *.o *.a *.oct *.$(EXT_DYN)

tar:
	tar czhf tools.tgz src/*.c include/*.h Makefile
