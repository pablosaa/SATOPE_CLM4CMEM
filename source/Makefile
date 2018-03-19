HOST = $(shell hostname)

ifeq "$(HOST)" "mysu.yolle"
export FORTRAN_COMPILER = gfortran
else
ifeq "$(HOST)" "juqueen3.fz-juelich.de"
export FORTRAN_COMPILER = xlf90
else
export FORTRAN_COMPILER = ifort
endif
endif

all:
	@echo "Compiling for " $(HOST)
ifdef $(STDALONE)
	+$(MAKE) -C cmem_v5.1/ STDALONE=1
else
	+$(MAKE) -C cmem_v5.1/
endif
