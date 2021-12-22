IDIR =../include
IDIRSYS = 
IDIRGNU = 
CC=gcc
CFLAGS=-I$(IDIR)
CXX=gcc
CCFLAGS=-I$(IDIR) -lm  -lfftw3 -lstdc++ -O3
EMCC=emcc
EMCCFLAGS=-I$(IDIR) -lfftw3 -L../fftw-3.3.10/.libs -s WASM=1 -s "EXPORTED_RUNTIME_METHODS=['cwrap']" -s EXPORTED_FUNCTIONS='["_mcLJ","_main"]' -O3 -Wpedantic -s ALLOW_MEMORY_GROWTH=1 -s ASSERTIONS=1 -s NO_EXIT_RUNTIME=1 

ODIR=obj
LDIR =../lib

_DEPS =
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

Linux=mclj mcljfft
Web=mclj.js mcljfft.js
linux: $(Linux)
web: $(Web)
all: $(Linux) $(Web)

OBJSHARED=atoms3d.o atomsbase.o erngls.o ernhs.o ernvdw.o mcboxbase.o mcnvtbase.o mcnvtfft.o prng.o timebase.o bar.o barhist.o ljs.o fftw.o grdR.o grdL.o
SHARED=atoms3d.cpp atomsbase.cpp erngls.cpp ernhs.cpp ernvdw.cpp mcboxbase.cpp mcnvtbase.cpp mcnvtfft.cpp prng.cpp timebase.cpp bar.c barhist.c ljs.c fftw.c grdR.c grdL.c

_OBJ = mclj.o $(OBJSHARED)
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_OBJFFT = mcljfft.o $(OBJSHARED)
OBJFFT = $(patsubst %,$(ODIR)/%,$(_OBJFFT))

$(ODIR)/%.o: %.c $(_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.cpp $(_DEPS)
	$(CXX) -c -o $@ $< $(CCFLAGS)

mclj: $(OBJ)
	$(CC) -o $@ $^ $(CCFLAGS)

mcljfft: $(OBJFFT)
	$(CC) -o $@ $^ $(CCFLAGS)

%.js: %.cpp $(SHARED) 
	$(EMCC) -o $@ $^ $(EMCCFLAGS)

.PHONY: clean deepclean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~

deepclean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~
	rm -f $(Linux) 
	rm -f $(Web) *.wasm
