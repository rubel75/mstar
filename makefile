# EXE OPTIONS:
#
# make
# make all (same as make)
# make clean
# make veryclean

# Edit to adjust for Fortran compiler and flags.

# Intel Fortran + MKL (debug)

#FC = ifort
#FCFLAGS = -I${MKLROOT}/include -g -traceback -check all -debug all -qopenmp -O0
#FLFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -qopenmp -ldl -g -traceback -check all -debug all -O0

# Intel Fortran + MKL (performance)

FC = ifort
FCFLAGS =  -I${MKLROOT}/include -qopenmp
FLFLAGS =  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -qopenmp -ldl

# GNU Fortran + MKL (performance)

#FC = gfortran
#FCFLAGS = -I${MKLROOT}/include -fopenmp -g -fbacktrace -ffpe-summary=none
#FLFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -fopenmp -ldl -g -fbacktrace -ffpe-summary=none

# GNU Fortran + OpenBLAS (performance)

#FC = gfortran
#OPENBLAS_ROOT = /you/path/to/openblas
#FCFLAGS = -I$(OPENBLAS_ROOT)/include -fopenmp -g -fbacktrace -ffpe-summary=none
#FLFLAGS = -L$(OPENBLAS_ROOT)/lib -lopenblas -lpthread -fopenmp -ldl -g -fbacktrace -ffpe-summary=none

# ~~~ Do not edit after that line ~~~

# Check math libraries variables
ifeq ($(strip $(MKLROOT)$(OPENBLAS_ROOT)),)
    $(error Neither MKLROOT (env) nor OPENBLAS_ROOT (internal) is set)
endif

PROGRAM = mstar

# source files and objects
SRCS = $(patsubst %.f90, %.o, $(wildcard *.f90)) \
    $(patsubst %.h, %.mod, $(wildcard *.h))


all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FLFLAGS) $(FLINK) -o $@ $^

%.o: %.f90
	$(FC) $(FCFLAGS) $(FOPT) -c $<

#%.mod: %.h
#	$(FC) $(FCFLAGS) -o $@ $<



# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
