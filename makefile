# EXE OPTIONS:
#
# make
# make all (same as make)
# make clean
# make veryclean

# Edit to adjust for Fortran compiler and flags.

# Intel Fortran

FC = ifort

#FCFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -qopenmp -ldl -g -traceback -check all -debug all
#FLFLAGS = -i8 -I${MKLROOT}/include -g -traceback -check all -debug all -qopenmp

FCFLAGS =  -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -qopenmp -ldl
FLFLAGS =  -i8 -I${MKLROOT}/include -qopenmp

# GNU Fortran

#FC = gfortran
#FCFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -fopenmp -ldl -g -fbacktrace -ffpe-summary=none
#FLFLAGS = -I${MKLROOT}/include -fopenmp -g -fbacktrace -ffpe-summary=none

# ~~~ Do not edit after that line ~~~

PROGRAM = mstar

# source files and objects
SRCS = $(patsubst %.f90, %.o, $(wildcard *.f90)) \
    $(patsubst %.h, %.mod, $(wildcard *.h))


all: $(PROGRAM)

$(PROGRAM): $(SRCS)
	$(FC) $(FCFLAGS) $(FLINK) -o $@ $^

%.o: %.f90
	$(FC) $(FLFLAGS) $(FOPT) -c $<

#%.mod: %.h
#	$(FC) $(FLFLAGS) -o $@ $<



# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAMS)
