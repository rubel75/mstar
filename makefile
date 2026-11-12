# EXE OPTIONS:
#
# make
# make all (same as make)
# make clean
# make veryclean
#
# Version information is generated automatically at build time and written
# to the output files and to the screen:
#
#     Git revision: YYYY-MM-DD-<12-character commit hash>
#
# If tracked files differ from HEAD, "-dirty" is appended. For a copy without
# Git metadata (e.g., a ZIP file downloaded from GitHub):
#
#     Git revision: unavailable; build date: YYYY-MM-DD

# Edit to adjust for Fortran compiler and flags.

# Intel Fortran + MKL (debug)

#FC = ifort
#FCFLAGS = -I${MKLROOT}/include -g -traceback -check all -debug all -O0
#FLFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -ldl -g -traceback -check all -debug all -O0

# Intel Fortran + MKL (performance)

FC = ifort
FCFLAGS =  -I${MKLROOT}/include -O2
FLFLAGS =  -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -ldl

# GNU Fortran + MKL (performance)

#FC = gfortran
#FCFLAGS = -I${MKLROOT}/include -O2 -g -fbacktrace -ffpe-summary=none
#FLFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -ldl -g -fbacktrace -ffpe-summary=none

# GNU Fortran + OpenBLAS (performance)

#FC = gfortran
#OPENBLAS_ROOT = /you/path/to/openblas
#FCFLAGS = -I$(OPENBLAS_ROOT)/include -O2 -g -fbacktrace -ffpe-summary=none
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

VERSION_SRC = version_mod.f90

# Modules must be compiled before the files that use them.
MODULE_OBJS = \
	precision_mod.o \
	version_mod.o \
	is_hermitian_mod.o \
	eigvz_mod.o \
	eigvs_mod.o \
	dgenen_mod.o \
	degenmass_mod.o \
	mommat_header_mod.o \
	read_waveder_header_mod.o \
	read_eigenvalues_vasp_mod.o \
	read_mommat_pij_vasp_mod.o \
	read_mommat_pij_wien2k_kpoint_mod.o

OBJS = $(MODULE_OBJS) mstar.o


all: $(PROGRAM)

# Libraries are listed after the objects that reference them.
$(PROGRAM): $(OBJS) | check-compiler
	$(FC) -o $@ $^ $(FLFLAGS) $(FLINK)

%.o: %.f90 | check-compiler
	$(FC) $(FCFLAGS) $(FOPT) -c $<

# Version information (as in BerryCPT)
#
# If this is a cloned Git repository, obtain the commit date and abbreviated
# commit hash from Git and append "-dirty" if tracked files differ from HEAD.
# Otherwise, record the build date. The file version_mod.f90 is replaced only
# when its content changes, to avoid unnecessary recompilation.

$(VERSION_SRC): FORCE
	@version=""; \
	if command -v git >/dev/null 2>&1 && \
	   [ -e .git ] && \
	   git rev-parse --is-inside-work-tree >/dev/null 2>&1; then \
		git_date=$$(git log -1 --format=%cs 2>/dev/null); \
		git_hash=$$(git rev-parse --short=12 HEAD 2>/dev/null); \
		if [ -n "$$git_date" ] && [ -n "$$git_hash" ]; then \
			dirty=""; \
			if ! git diff --quiet HEAD -- 2>/dev/null; then \
				dirty="-dirty"; \
			fi; \
			version="Git revision: $$git_date-$$git_hash$$dirty"; \
		fi; \
	fi; \
	if [ -z "$$version" ]; then \
		build_date=$$(date +%Y-%m-%d); \
		version="Git revision: unavailable; build date: $$build_date"; \
	fi; \
	{ \
		echo 'module version_mod'; \
		echo '    implicit none'; \
		echo '    private'; \
		echo '    public :: version_string'; \
		echo ''; \
		printf '    character(len=*), parameter :: version_string = "%s"\n' "$$version"; \
		echo ''; \
		echo 'end module version_mod'; \
	} > $(VERSION_SRC).tmp; \
	if [ ! -f $(VERSION_SRC) ] || \
	   ! cmp -s $(VERSION_SRC).tmp $(VERSION_SRC); then \
		mv $(VERSION_SRC).tmp $(VERSION_SRC); \
	else \
		rm -f $(VERSION_SRC).tmp; \
	fi

# Module dependencies
is_hermitian_mod.o: precision_mod.o
eigvz_mod.o: precision_mod.o
eigvs_mod.o: precision_mod.o
dgenen_mod.o: precision_mod.o
degenmass_mod.o: precision_mod.o eigvz_mod.o is_hermitian_mod.o
read_waveder_header_mod.o: precision_mod.o
read_eigenvalues_vasp_mod.o: precision_mod.o
read_mommat_pij_vasp_mod.o: precision_mod.o
read_mommat_pij_wien2k_kpoint_mod.o: precision_mod.o
mstar.o: $(MODULE_OBJS)


# Utility targets
.PHONY: all clean veryclean check-compiler FORCE

FORCE:

check-compiler:
	@if ! command -v "$(firstword $(FC))" >/dev/null 2>&1; then \
		echo "ERROR: Fortran compiler '$(firstword $(FC))' was not found."; \
		echo "       Set FC in the makefile to an available compiler."; \
		exit 1; \
	fi

clean:
	rm -f *.o *.mod *.MOD *__genmod.*

veryclean: clean
	rm -f *~ $(PROGRAM) $(VERSION_SRC) $(VERSION_SRC).tmp
