# Makefile for wrappers
# COSMOPMC/MonteCarlo/src/wrappers/Makefile
# Martin Kilbinger, Karim Benabed 2007-2009

include $(COSMOPMC)/Makefile.main


CFLAGS      = -Iinclude -I$(GSL)/include -I$(FFTW)/include -I$(TOOLS)/include \
	      $(IINCDIRS) $(INC_PMC) $(OPTIM) $(WARNINGS) \
              $(ARCHFLAGS) $(FPIC) $(DEBUG) $(MKCFLAGS)
LFLAGS      = -L$(TOOLS) $(LLIBDIRS) $(LLIBS) -lm -lmvdens -lerrorio -lpar -lconf $(FPIC)
LIBS        = $(LIBERROR) libmvdens.$(EXT)

# Include all files for RELEASE version, and only selected files for DEVEL version
ifeq ($(STATUS), "RELEASED")
   OBJ         := $(patsubst src/%.c,%.o,$(wildcard src/*.c))
endif
ifeq ($(STATUS), "DEVEL")
  OBJ          = bao.o halo.o sn.o lens.o timexec.o wmap.o xcorr.o nz.o wrappers.o init_wrappers.o
endif

# New: the following lines were only in the 'RELEASE' loop, for an unknown reason
ifeq ($(DOWMAP),1)
   CFLAGS  += -DDOWMAP
endif


vpath %.h ../wrappers/include
vpath %.c ../wrappers/src



all: $(OBJ) libparam.a libwr.$(EXT) libpar.$(EXT)


%.o: %.c
	$(CC) $< -c $(CFLAGS)

libparam.$(EXT_DYN):
	$(CC) $(SHARED) -o $ $@ param.o $(LFLAGS)

libparam.a: param.o
	$(AR) -r $@ $^

libwr.$(EXT_DYN): $(OBJ) $(LIBS)
	$(CC) $(SHARED) -o $ $@ $(OBJ) $(LFLAGS)

libwr.a: $(OBJ)
	$(AR) -r $@ $^

$(LIBERROR) libmvdens.$(EXT) :
ifeq ($(PMCLIB_EXT),0)
	cd $(TOOLS) && $(MAKE) $@
endif

libpar.$(EXT) libconf.$(EXT) :
	cd $(TOOLS) && $(MAKE) $@


libcm.$(EXT):
	cd $(MF)/obj && $(MAKE)	

%.h:


PHONY: clean tar

clean:
	rm -rf *.o *~ lib*.$(EXT_DYN) libwr*.a libparam.a libpar.a

tar:
	tar czf wrapper.tgz src/*.c include/*.h Makefile
