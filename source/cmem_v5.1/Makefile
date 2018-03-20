# Makefile for CMEM and CLM SATELLITE OPERATOR sources
#

# Options to compile are as "cmem" or as a example "TESTBED"
# WHERE:
# SOURCE = cmem_main.F90 is for the stand-alone version
# SOURCE = TESTBED.F90 is for the subroutine version, "TESTBED.F90" can
#          be any code which need to get the TBs from CMEM SAT OPERATOR
#


ifdef MUTTER
	SOURCE = $(MUTTER)
	TMP    = $(MUTTER:.F90=)
	EXECFI = $(notdir $(TMP))
	# EXECFI = TESTBED 
	# SOURCE = TESTBED.F90 
else
	EXECFI = cmem
	SOURCE = cmem_main.F90
endif

OBJECT = $(notdir $(SOURCE:.F90=.o))

BIN_PATH = ../../bin
NETCDF_INCLUDE_PATHS = -I/usr/include 
NETCDF_LIBRARY_PATHS = -L/usr/lib64 
NETCDF_LIBRARY_NAMES = -lnetcdff
CLM4CMEM_PATH = ../clm4cmem/

RM = rm -f
# PSG adding: the move command so executable is in ../bin
MV = mv $(EXECFI) $(BIN_PATH)

$(EXECFI) : $(OBJECT)
	${FORTRAN_COMPILER} -o $(EXECFI) parkind1.o yomlun_ifsaux.o yomlun.o abor1.o \
yomcmemnetcdf.o yomcmempar.o \
io_cmemnetcdf.o yomcmemsoil.o \
yomcmemveg.o yomcmematm.o yomcmemfields.o cmem_setup.o  vegtable.o \
atm_sub.o \
rdcmemnetcdfinfo.o rdcmemnetcdf.o rdcmemasciiinfo.o rdcmemascii.o cmem_init.o dielice_sub.o \
dielwat_sub.o dielsoil_sub.o teff_sub.o fresnel_sub.o wilheit_sub.o \
rghref_sub.o veg_sub.o cmem_soil.o cmem_veg.o cmem_snow.o cmem_atm.o \
cmem_rtm.o wrcmemascii.o wrcmemnetcdf.o \
rdclm_wrcmem.o toolbox_clm4cmem.o toolbox_SatOperator.o \
SAT_CLM4CMEM.o $(OBJECT) \
${NETCDF_LIBRARY_PATHS} ${NETCDF_LIBRARY_NAMES}


	$(RM) *o
	$(RM) *mod
	$(MV)

clean:
	$(RM) *o
	$(RM) *mod

#toolbox_clm4cmem.o : $(CLM4CMEM_PATH)toolbox_clm4cmem.F90
#	${FORTRAN_COMPILER} -c $(CLM4CMEM_PATH)toolbox_clm4cmem.F90

$(OBJECT) : $(SOURCE)
	${FORTRAN_COMPILER} -c ${NETCDF_INCLUDE_PATHS} ${GRIBAPI_INCLUDE_PATHS} parkind1.F90 yomlun_ifsaux.F90 yomlun.F90 abor1.F90 \
	yomcmemnetcdf.F90 yomcmempar.F90 \
	io_cmemnetcdf.F90 yomcmemsoil.F90 \
	yomcmemveg.F90 yomcmematm.F90 yomcmemfields.F90 \
	vegtable.F90 atm_sub.F90 \
	$(CLM4CMEM_PATH)toolbox_clm4cmem.F90 \
	$(CLM4CMEM_PATH)rdclm_wrcmem.F90 \
	cmem_setup.F90 rdcmemnetcdfinfo.F90 rdcmemnetcdf.F90 rdcmemasciiinfo.F90 \
	rdcmemascii.F90 cmem_init.F90 dielice_sub.F90 \
	dielwat_sub.F90 dielsoil_sub.F90 teff_sub.F90 \
	fresnel_sub.F90 wilheit_sub.F90 \
	rghref_sub.F90 veg_sub.F90 cmem_soil.F90 cmem_veg.F90 \
	cmem_snow.F90 cmem_atm.F90 \
	cmem_rtm.F90 wrcmemascii.F90 wrcmemnetcdf.F90 \
	$(CLM4CMEM_PATH)toolbox_SatOperator.F90 \
	$(CLM4CMEM_PATH)SAT_CLM4CMEM.F90 $(SOURCE)


# Non used variables:
#GRIBAPI_INCLUDE_PATHS = -I/usr/local/apps/grib_api/1.12.3/GNU/4.8.1/include
#GRIBAPI_LIBRARY_PATHS = -L/usr/local/apps/grib_api/1.12.3/GNU/4.8.1/lib -L/usr/local/apps/grib_api/1.12.3/GNU/4.8.1/lib
#GRIBAPI_LIBRARY_NAMES = -lgrib_api_f90 -lgrib_api -ljasper

# PSG extracting: yomcmemgribapi.F90, io_cmemgribapi.F90, rdcmemgribapi.F90,  rdcmemgribapiinfo.F90, wrcmemgribapi.F90
