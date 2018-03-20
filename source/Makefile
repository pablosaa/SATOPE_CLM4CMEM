HOST = $(shell hostname)

ifeq "$(HOST)" "mysu.yolle"
export FORTRAN_COMPILER = gfortran
else
ifeq "$(HOST)" "juqueen.fz-juelich.de"
export FORTRAN_COMPILER = xlf90
else
export FORTRAN_COMPILER = ifort
endif
endif

all:
	@echo "Compiling for " $(HOST)
ifdef MUTTER
	# Fot Subroutine CMEM+SAT_OPE called from main program MUTTER
	+$(MAKE) -C cmem_v5.1/ MUTTER=$(MUTTER)
	@echo "*** Satellite Operator has been created to be called from" $(MUTTER)
else
	# For Stand-alone CMEM+SAT_OPE compilation
	+$(MAKE) -C cmem_v5.1/
	@echo "*** CMEM stand-alone has been created!"
endif

clean:
	+$(MAKE) -C cmem_v5.1/ clean


help:
	@echo "USAGE:"
	@echo "To generates the CMEM stand-alone version, execute:"
	@echo " > make "
	@echo "or..."
	@echo "To generates the SAT OPERATOR subroutine version which will be "
	@echo "called from any program (in this example TESTBED.F90), execute:"
	@echo " > make MUTTER=/full_path/to_mother/code/TESTBED.F90"
	@echo "note that MUTTER must be indicated with the full path of the example code TESTBED.F90".
	@echo "Alternativelly if the mother source code e.g. TESTBED.F90 is located in actual dir: "
	@echo " > make MUTTER=$$´{PWD}/TESTPED.F90"
	@echo "might works as well..."
	@echo "Sometimes it is better to clean object and modules, before new compilation:"
	@echo " > make clen "
	@echo ":)"


