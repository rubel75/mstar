# EXE OPTIONS:
#
# make
# make all (same as make)
# make clean
# make veryclean

# Edit to adjust for Fortran compiler and flags.

# Intel Fortran + MKL (debug)

#FC = ifort
#FCFLAGS = -I${MKLROOT}/include -g -traceback -check all -debug all -O0
#FLFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -ldl -g -traceback -check all -debug all -O0

# Intel Fortran + MKL (performance)

FC = ifort
FCFLAGS =  -I${MKLROOT}/include
FLFLAGS =  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -ldl

# GNU Fortran + MKL (performance)

#FC = gfortran
#FCFLAGS = -I${MKLROOT}/include -g -fbacktrace -ffpe-summary=none
#FLFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -ldl -g -fbacktrace -ffpe-summary=none

# GNU Fortran + OpenBLAS (performance)

#FC = gfortran
#OPENBLAS_ROOT = /you/path/to/openblas
#FCFLAGS = -I$(OPENBLAS_ROOT)/include -g -fbacktrace -ffpe-summary=none
#FLFLAGS = -L$(OPENBLAS_ROOT)/lib -lopenblas -lpthread -ldl -g -fbacktrace -ffpe-summary=none

# ~~~ Do not edit after that line ~~~

# Check math libraries variables
# (not required for the clean and veryclean targets)
ifeq ($(filter clean veryclean,$(MAKECMDGOALS)),)
ifeq ($(strip $(MKLROOT)$(OPENBLAS_ROOT)),)
    $(error Neither MKLROOT (env) nor OPENBLAS_ROOT (internal) is set)
endif
endif

PROGRAM = mstar

# Modules must be compiled before the files that use them.
MODULE_OBJS = \
	precision_mod.o \
	is_hermitian_mod.o \
	eigvz_mod.o \
	eigvs_mod.o \
	dgenen_mod.o \
	degenmass_mod.o \
	read_waveder_header_mod.o \
	read_eigenvalues_vasp_mod.o \
	read_mommat_pij_vasp_mod.o

OTHER_OBJS = \
	kdelta.o \
	read_mommat_nb.o \
	read_mommat_pij.o \
	read_numlines.o

OBJS = $(MODULE_OBJS) $(OTHER_OBJS) mstar.o


all: $(PROGRAM)

# Libraries are listed after the objects that reference them.
$(PROGRAM): $(OBJS)
	$(FC) -o $@ $^ $(FLFLAGS) $(FLINK)

%.o: %.f90
	$(FC) $(FCFLAGS) $(FOPT) -c $<

# Module dependencies
is_hermitian_mod.o: precision_mod.o
eigvz_mod.o: precision_mod.o
eigvs_mod.o: precision_mod.o
dgenen_mod.o: precision_mod.o
degenmass_mod.o: precision_mod.o eigvz_mod.o is_hermitian_mod.o
read_waveder_header_mod.o: precision_mod.o
read_eigenvalues_vasp_mod.o: precision_mod.o
read_mommat_pij_vasp_mod.o: precision_mod.o
mstar.o: $(MODULE_OBJS)


# Utility targets
.PHONY: all clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAM)
