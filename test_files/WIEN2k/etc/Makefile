# make all			-- run whole sequence starting with cleaning
# make clean		-- clean data
# make init			-- initialize calculation
# make run			-- run calculation
# make postproc		-- extract results
#
# Files needed: Bi2Se3.struct
# Python is needed (module load python)

# speed of light in at. units: 137 (regular SOC) or 300 (reduced SOC)
CLIGHT=137
BERRYCPTROOT=~/berrycpt
CASE=$(notdir $(shell pwd))


all: clean init run postproc

init: CoSi.struct | check-env
	date > date
	@echo "Initializing calculation"
	@echo CASE=$(CASE)
	cp CoSi.struct $(CASE).struct
	x sgroup
	cp $(CASE).struct_sgroup $(CASE).struct
	init_lapw -b -rkmax 8 -vxc 13 -red 3 -prec 2n -nohdlo
	export EDITOR=cat; \
	echo -e "0 0 1\n\n\n\nN" | init_so_lapw; \
	echo -e "-1\n0.15\n1" | x kgen -so
	touch $@ # create dummy target file

run: init CoSi.klist_band | check-env
	date >> date
	ln -s -f $(WIENROOT)/SRC_lapwso/lapwso-c$(CLIGHT) $(WIENROOT)/lapwso # reduced SOC
	ls -l $(WIENROOT)/lapwso # check SOC strength
	export OMP_NUM_THREADS=$(SLURM_CPUS_PER_TASK); \
	run_lapw -ec 0.0001 -cc 0.001 -so
	save_lapw -d SOC # save SCF
#	cp SOC/$(CASE).in1c $(CASE).in1
	sed -i 's/4   -9\.0       5\.0/4   -9.0      10.0/g' $(CASE).in1c # Emax=10 Ry for LAPW1
	sed -i 's/^-10  1\.9/-10  7.0/g' $(CASE).inso # Emax=7 Ry for LAPWSO
	cp CoSi.klist_band $(CASE).klist_band
	export OMP_NUM_THREADS=$(SLURM_CPUS_PER_TASK); \
	x lapw1 -band
	export OMP_NUM_THREADS=$(SLURM_CPUS_PER_TASK); \
	x lapwso
	ln -s -f $(WIENROOT)/SRC_lapwso/lapwso-c137 $(WIENROOT)/lapwso # set back SOC to normal
	touch $@ # create dummy target file

postproc: run | check-env check-berrycpt
	export OMP_NUM_THREADS=$(SLURM_CPUS_PER_TASK); \
	rm $(CASE).vspup
	ln -s $(CASE).vsp $(CASE).vspup
	rm $(CASE).vspdn
	ln -s $(CASE).vsp $(CASE).vspdn
	ln -s $(CASE).vectorso $(CASE).vectorsoup
	cp $(WIENROOT)/SRC_templates/case.inop $(CASE).inop # get input file template
	sed -i 's/^OFF/ON/g' $(CASE).inop # enable writing of momentum matrix elements in *.mommat2*
	sed -i 's/^-5\.0 3\.0/-5.0 7.0/g' $(CASE).inop # extend Emax = 7 Ry
	export OMP_NUM_THREADS=$(SLURM_CPUS_PER_TASK); \
	x optic -so -up
	export OMP_NUM_THREADS=$(SLURM_CPUS_PER_TASK); \
	$(BERRYCPTROOT)/berrycpt $(CASE).mommat2 -nvb 78 -so
	touch $@ # create dummy target file

clean:
	rm -rf init run postproc # dummy targets
	rm -rf :* *.day* $(CASE).* *.def *.error *clm* SOC *.dat *.pdf *.csv *.png kvectors date
	rm -rf .[!.]* ..?* # files .xx

.PHONY: all clean check-env

# check environment variables
check-env:
ifndef WIENROOT
	$(error WIENROOT environment variable is undefined)
endif

# check berrycpt code
check-berrycpt:
ifeq ("$(wildcard $(BERRYCPTROOT)/berrycpt)","")
    $(error Unable to resolve the path to berrycpt: $(BERRYCPTROOT)/berrycpt)
endif
