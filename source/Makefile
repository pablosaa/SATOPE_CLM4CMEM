HOST = $(shell hostname)

ifeq ("$(HOST)", $(filter "$(HOST)", "mysu.yolle" "r233"))
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
	@echo "otherwise..."
	@echo "To generates the SAT OPERATOR subroutine version which will be "
	@echo "called from any other program (in this example TESTBED.F90), execute:"
	@echo " > make MUTTER=/full_path/to_mother/code/TESTBED.F90"
	@echo "note that MUTTER must be indicated with the full path of the example code TESTBED.F90".
	@echo "Alternativelly if the mother source code e.g. TESTBED.F90 is located in the actual directory: "
	@echo " > make MUTTER=$$´{PWD}/TESTPED.F90"
	@echo "that might works as well..."
	@echo "Sometimes it is better to clean object and module files before new compilation,"
	@echo "in that case run:"
	@echo " > make clean "
	@echo ":)"


