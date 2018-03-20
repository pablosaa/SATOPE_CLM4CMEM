# SATOPE_CLM4CMEM
Code for a satellite simulator alike SMOS or SMAP by means of CMEM as radiative transfer and CLM as input for land surface model.

## STRUCTURE ##
	./SATOPE_CLM4CMEM/
		|______ bin/
		|______ forcing/
		|______ inputlst/
		|______ output/
		|______ source/
			|_____ cmem_v5.1/
			|_____ clm4cmem/

Where:
*/bin* directory is where executables are placed

*/forcing* is the directory where forcing files must be located, e.g. Satellite orbits, classical CMEM NetCDF input files (optional, not necessary) , CLM input files (optional)

*/inputlst* this directory host the files with parameters defined as ASCII "input" following the classical CMEM "input", with the difference that many can be storaged with different names and specified them respectively for specific runnings.

*/output* Output NetCDF files (either for classical CMEM or for the satellite operator) will be located in this directory after running the code

*/source* this directory is separated in two parts: */cmem__v5.1* where the main CMEM source code is located (this however contains part of source codes which have been modified from the official CMEM distribution), and */clm4cmem* where all codes for the tools necessary to apply the satellite operator are located.

## CONFIGURE AND BUILDING ##
The most important feature to configure is the ``FORTRAN_COMPILER`` to use, this can be done by tuning the environment variable ``FORTRAN_COMPILER`` in the Makefile located at /source.

Additionally, the location for the NetCDF ``include`` and ``lib`` directories might be needed to be tuned, that can be done in the Makefile located at /source/cmem__v5.1

The building process offers two options, namely as Stand-alone executable or as Subroutive to be called from other codes. 
Go to the source directory and see help for details, e.g. type: ``> make help``

1. To build the Stand-alone version, type: ``> make``
	this will produce the executable cmem at the bin directory.
2. To build the Subroutine version, which will be used under a generic (as example) code file ``TESTBED.F90``, type: ``> make MUTTER=/path/where/codes-is/TESTBED.F90``
this will create the executable ``TESTBED`` at the */bin* directory, this executable includes the Satelletite Operator as subroutine within.

## RUNNING ##
The abobe mentioned two flavours to run the Satellite Operator, i.e. Stand-alone version or Subroutine version can be running as follow:
1. Stand-alone: executable file "cmem"

		USAGE: > ./cmem CLM_input_file.nc [input_parameters_file]
	Where [input_parameters_file] is optional and indicates the typical CMEM ASCII input file defining the CMEM parameters. By default it will search the file ``input`` in the */inputlst* directory.

2. Subroutine version: This is for integration of the satellite operator with any software who needs to call the operator. The name of subroutine is ``sat__clm4cmem()`` 

		USAGE: call sat_clm4cmem(CLM_input_file,input_parameters,SAT)
where ``CLM_input_file`` must be a variable type character containing the full-path of the CLM input file to work.
		``input_parameters`` (optional) must be a variable type character containing the name of the ASCII input file with the CMEM parameters to be used, by default the code looks for a file named ``input`` in */inputlst* directory.
	
_EXAMPLE_: whithin some fortran code the Satellite Operator can be invocated as illustrated in the following example:

	PROGRAM TESTBED

	use clm4cmem, only: SATELLITE
	implicit none

	character(len=300) :: CLM_input_file, input_parameters
	type(SATELLITE) :: SAT

	CLM_input_file = "/run/data/3022-C347/clmoas.clm2.h0.2008-03-15-00900.nc"
	input_parameters = "input_test1"
	call sat_clm4cmem(trim(CLM_input_file),trim(input_parameters),SAT)

	END PROGRAM

In the example above, the program ``TESTBED`` will call the Satellite Simulator passing the CLM data file */run/data/3022-C347/clmoas.clm2.h0.2008-03-15-00900.nc* as input to be precessed, with CMEM parameters specified by the ASCII file name *input_test1* which must be located at */inputlst* directory.

The variable structure ``SAT`` contains the output of the satellite operator as well as the satellite information provided as input (see below in parameters). This variable ``SAT`` can be used for any purposes by other procedures.

The variable type ``SATELLITE`` is a structure defined at */source/clm4cmem/toolbox_clm4cmem.F90* and it is comprised by the following fields:

	type SATELLITE
		character(len=10) :: name
		real(kind=JPRM) :: orbit
		real(kind=JPRM) :: antenna, azimuth, wavelength
		real(kind=JPRM), allocatable, dimension(:) :: theta, incl_foprt, time
		real(kind=JPRM), allocatable, dimension(:) :: lon_foprt, lat_foprt
		real(kind=JPRM), allocatable, dimension(:,:,:,:) :: TBSAT_HV
		character(len=300) :: OrbitFileName
	end type SATELLITE
Where:
* name: string with Sensor's name for referece e.g. "SMOS"
* orbit: Orbit altitude of Satellite [km]
* antenna: Sensor's antenna diameter [m]
* wavelength: Sensor's wavelength [m]
* theta: Vector with incidence angles [deg] e.g (/20, 30, 40, 50/)
* lon_foprt, lat_foprt: longitude and latitude of footprints [deg]
* azimuth: Along-track Satellite's orbit azimuth [deg] from North ClockWise
* incl_foprt: Footprint inclination angle [deg] from Equator CounterClockWise (E-CCW) 
* TBSAT_HV: Bightness Temperature [K] (pixel,polarization,incidence_angle,time)
* OrbitFileName: string with file name of the NetCDF containing the Satellite information (i.e. footprint coordinates Longitude, Latitude, Inclination, wavelength, etc.).

## CONFIGURATION PARAMETERS ##

Parameter configuration is done by the usual ASCII input file, but now located at inputlst directory. The configuration files include now new parameters like:

* ``CFINOUT='clm'``, for specifying CLM as processing input file
* ``JPHISTLEV=4``,   to produce NetCDF level-4 output which is only Satellite Simulator data. This can also be from 1 to 3 (as typical CMEM) but also can be 5 for outputs as level 1 to 4, or 6 for no NetCDF output and the information is keeped only in memory only.
* ``INPUTSATINFO='../forcing/SMOS_L1orbit_neckar.nc'``, indicating the Satellite information to be used as NetCDF format.

...
 

