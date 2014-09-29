# compile executables each from a single source

FC       = gfortran
FFLAGS   = -ffree-line-length-none 

SOURCES  := $(wildcard *.f90)
BINARIES := $(SOURCES:.f90=.exe)

all: $(BINARIES)

%.exe: %.o
	$(FC) $(FFLAGS) $< -o $@

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f *.exe *.o

