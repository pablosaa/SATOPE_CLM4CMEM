# SATOPE_CLM4CMEM
Code for a satellite simulator alike SMOS or SMAP by means of CMEM as radiative transfer and CLM as input for land surface model.

Published at Zenodo under: 
[![DOI number](https://zenodo.org/badge/DOI/10.5281/zenodo.3516115.svg)](https://zenodo.org/record/3516115)

## STRUCTURE ##
	./SATOPE_CLM4CMEM/
		|______ bin/
		|______ forcing/
		|______ inputlst/
		|______ output/
		|______ scripts/
		|______ source/
			|_____ cmem_v5.1/
			|_____ clm4cmem/

Where:
*/bin* directory is where executable binaries are placed

*/forcing* is the directory where forcing files must be located, e.g. Satellite orbits, classical CMEM NetCDF input files (optional, not necessary) , CLM input files (optional)

*/inputlst* this directory host the files with parameters defined as ASCII "input" following the classical CMEM "input", with the difference that many can be created with different names and holding different configurations, then one of them can be specified as input parameter file for running CMEM with that specific parametrization.

*/output* Output NetCDF files (either for classical CMEM or for the satellite operator) will be located in this directory after running the code

*/scripts* Directory containing some useful MATLAB/GNU_Octave scripts to deal with the data. e.g. plotting multi-angle multi-files Level 4 CMEM output NetCDF data, plotting Level 1 CMEM output NetCDF data, scripts to create the Satellite input NetCDF information (see below), etc.

*/source* this directory is separated in two parts: */cmem_v5.1* where the main CMEM source code is located (this however contains part of source codes which have been modified from the official CMEM distribution), and */clm4cmem* where all codes for the tools necessary to apply the satellite operator are located.

## CONFIGURE AND BUILDING EXECUTABLE ##
The most important feature to configure is the environment variable ``FORTRAN_COMPILER`` that defines which compiler to use, this can be done by tuning ``FORTRAN_COMPILER`` in the *Makefile* located at */source*.

Additionally, the location for the NetCDF ``include`` and ``lib`` directories might be needed to be adjusted (system dependent), which can be done editing the ``Makefile`` located at */source/cmem_v5.1*

Compilation and running have been tested in Linux machines with compilers like **gfrotran** (Linux station and CLUMA clusted at Meteorological Institute) and **ifortran** (JURECA supercomputer at Jülich Computing Center).

The building process offers the possibility to generate the CMEM satellite operator executable in two flavors, namely as Stand-alone executable or as Subroutine to be called from other program.

Go to the source directory and see help for details, e.g. in command line type: ``> make help``

1. To build the Stand-alone version, type: ``> make``
	this will produce the executable cmem at the bin directory.
2. To build the Subroutine version, which can be used by any generic code. Here as example a code named ``TESTBED.F90`` is used, type: ``> make MUTTER=/path/where/codes-is/TESTBED.F90``
this will create the executable ``TESTBED`` at the */bin* directory, this executable includes the Satellite Operator as subroutine within.

## RUNNING ##
The above mentioned two flavours to run the Satellite Operator, i.e. Stand-alone version or Subroutine version can be running as follow:
1. Stand-alone: executable file "cmem"

		USAGE: > ./cmem CLM_input_file.nc [input_parameters_file]
	Where [input_parameters_file] is optional and indicates the typical CMEM ASCII input file defining the CMEM parameters. By default it will search the file ``input`` in the */inputlst* directory.

_EXAMPLE_: from command line go to the SATOPE_CLM4CMEM/bin directory and run the program as simple as:

		> ./cmem /tmp/clm/clmmoas/3022-C347/clmoas.clm2.h0.2008-03-15-00900.nc

this will run CMEM considering as forcing CLM NetCDF file ``clmoas.clm2.h0.2008-03-15-00900.nc`` located at the directory ``/tmp/clm/clmmoas/3022-C347/``. The ancillary NetCDF surface-file corresponding to the CLM file must also be located at that directory (optionally other location can be specified in the code).

2. Subroutine version: This is for integration of the satellite operator with any software who needs to call the operator. The name of subroutine is ``sat_clm4cmem()``

		USAGE: call sat_clm4cmem(CLM_input_file,input_parameters,SAT)
where ``CLM_input_file`` must be a variable type character containing the full-path of the NetCDF CLM forcing file to process.
		``input_parameters`` (optional) must be a variable type character containing the name of the ASCII input file with the CMEM parameters to be used, by default the code looks for a file named ``input`` in */inputlst* directory.

_EXAMPLE_: within any Fortran or C code the Satellite Operator can be accessed as illustrated in the following example:

		PROGRAM TESTBED

		use clm4cmem, only: SATELLITE
		implicit none

		character(len=300) :: CLM_input_file, input_parameters
		type(SATELLITE) :: SAT

		CLM_input_file = "/run/data/3022-C347/clmoas.clm2.h0.2008-03-15-00900.nc"
		input_parameters = "input_test1"
		call sat_clm4cmem(trim(CLM_input_file),trim(input_parameters),SAT)

		! form here on do whatever is needed with SAT type data...

		END PROGRAM

In the example above, the program ``TESTBED`` will call the Satellite Simulator passing the CLM data file */run/data/3022-C347/clmoas.clm2.h0.2008-03-15-00900.nc* as input to be processed, with CMEM parameters specified by the ASCII file name *input_test1* which must be located at */inputlst* directory.

The variable structure ``SAT`` (declared with ``type(SATELLITE) :: SAT``) contains the output of the satellite operator as well as the satellite information provided as sensor information (see below in [Configuration Parameters](#cmem-configuration-parameters)). This variable ``SAT`` can be used for any purposes by other procedures.

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
* name: string with Sensor's name for reference e.g. "SMOS"
* orbit: Orbit altitude of Satellite [km]
* antenna: Sensor's antenna diameter [m]
* wavelength: Sensor's wavelength [m]
* theta: Vector with incidence angles [deg] e.g (/20, 30, 40, 50/)
* lon_foprt, lat_foprt: longitude and latitude of footprints [deg]
* azimuth: Along-track Satellite's orbit azimuth [deg] from North ClockWise
* incl_foprt: Footprint inclination angle [deg] from Equator CounterClockWise (E-CCW)
* TBSAT_HV: Brightness Temperature [K] (pixel,polarization,incidence_angle,time)
* OrbitFileName: string with file name of the NetCDF containing the Satellite information (i.e. footprint coordinates Longitude, Latitude, Inclination, wavelength, etc.).

### Satellite Orbit Information ###

In order for the satellite operator to have information about the characteristics of the satellite's sensor and orbit, a NetCDF file containing basic information is needed. This NetCDF file is assigned to the field name ``OrbitFileName`` in the Fortran ``SATELLITE`` type variable (see section above). The name of the NetCDF file containing the desired satellite information is then passed to CMEM by  the configuration parameter ``INPUTSATINFO`` in the ASCII input file (see section below).

A typical NetCDf satellite information file, typically located at ``SATOPE_CLM4CMEM/forcing/`` directory, must have the following minimum structure:
		> ncdump -h SATOPE_CLM4CMEM/forcing/SMOS_L1orbit_neckar.nc

		netcdf SMOS_L1orbit_neckar {
		dimensions:
			NINC = 3 ;
			NPIXEL = 320 ;
		variables:
			float THETA_INC(NINC) ;
				THETA_INC:unit = "degrees_nadir" ;
				THETA_INC:long_name = "Incidence Angle" ;
				THETA_INC:short_name = "theta_inc" ;
			float LONG(NPIXEL) ;
				LONG:unit = "degree_east" ;
				LONG:long_name = "FootPrint Longitude" ;
				LONG:short_name = "foprt_lon" ;
			float LATI(NPIXEL) ;
				LATI:unit = "degree_north" ;
				LATI:long_name = "FootPrint Latitude" ;
				LATI:short_name = "foprt_lat" ;
			float INCLI(NPIXEL) ;
				INCLI:unit = "degree_equator" ;
				INCLI:long_name = "FootPrint Inclination" ;
				INCLI:short_name = "foprt_incl" ;

		// global attributes:
			:SATELLITE_name = "SMOS" ;
			:Orbit_altitude_km = 758. ;
			:Orbit_azimuth_deg = 98.1 ;
			:SENSOR_antenna_m = 7.5 ;
			:SENSOR_wavelength_m = 0.21 ;
			:creator contact = "username@uni-bonn.de" ;
		}

#### Dimensions and Variables Description ####

The NetCDF variables have two dimensions, ``NINC`` for the number of incidence angles and ``NPIXEL`` for the number of footprints to consider. In the example above the file has three incidence angles stored at the variable ``THETA_INC`` e.g. [30 40 50] in Degrees.

The variables with the ``NPIXEL`` dimension are ``LONG``, ``LATI`` and ``INCLI`` for the Longitude, Latitude and Footprint's inclination angle respectively. The inclination angle is given in degrees counter-clockwise from the Equator, while Longitude (Latitude) is in degrees positive for East (North).

The following Global attributes ``Orbit_altitude_km``, ``Orbit_azimuth_deg``, ``SENSOR_antenna_m`` and ``SENSOR_wavelength_m`` are very important for the satellite operator, and are sensor/satellite dependent. ``Orbit_altitude_km`` is a nominal satellite altitude for the region of study, ``Orbit_azimuth_deg`` is the orbit inclination translated to azimuth angle (i.e. degrees clockwise from North). ``SENSOR_antenna_m`` is the sensor's antenna nominal diameter in meters (i.e. for SMAP radiometer 6 m), and ``SENSOR_wavelength_m`` is the sensor's wavelength in meters (e.g. for typical L-band 0.21 m).

Additionally other attributes can be given as a reference for the file, e.g. ``SATELLITE_name`` as a name of the satellite/sensor and ``creator contact`` as a reference for the creator of the file or project.



## CMEM CONFIGURATION PARAMETERS ##

Parameter configuration is done by the usual ASCII input file, but now located at */inputlst* directory. The configuration files include now some new parameters and/or options:

* ``CFINOUT='clm'``, new option for specifying CLM as processing input file,
* ``JPHISTLEV=4``,   new option to produce NetCDF level-4 output which is only Satellite Simulator data. This can also be from 1 to 3 (as typical CMEM) but also can be 5 for outputs as level 1 to 4, or 6 for no output and the CMEM output information is kept only in memory only (which can be used by third-party codes),
* ``INPUTSATINFO='../forcing/SMOS_L1orbit_neckar.nc'``, new parameter indicating the Satellite information to be used as NetCDF format,
* ``INDEXTIME= 30,2,36`` -> new parameter to indicate which time indexes from the CLM input data to take into account, the format follows a sequence of three integers comprising (initial-index, total-number, span). In the given example (30,2,36) it means two time points with the first time having the index 30 and the second the index 30+36.

## SCRIPTS ##
In the directory ``SATOPE_CLM4CMEM/scripts`` few MATLAB/GNU Octave friendly scripts are located to help with the edition and visualization of basic input/output simulation related data for the CMEM-Satellite-Operator scheme. The following are some of them:

### Edit Satellite Info Input GUI ###
``edit_info_SAT4CMEM.m`` This is a simple GUI to manage the parameters contained in a NetCDF satellite information file like the ``SMOS_L1orbit_neckar.nc`` which is for SMOS. The GUI allows to load an existing NetCDF file, edit it or change the parameters and save it into a new file.

When run for the first time, the script only shows some values as example. To start working, it is suggested to load an existing file from ``SATOPE_CLM4CMEM/forcing`` directory.

The GUI is shown in Figure xxx where two main frames can be observed: the first corresponding to the global attributes and the second frame for the variables. The global attributes can be edited/changed directly, but the variables (e.g. ``THETA_INC``, ``LONG``, ``LATI`` and ``INCLI``) can be either edited directly or just indicated by the name of variable containing the data in the MATLAB/GNU Octave workspace.

_to do: include picture of the GUI_

#### SAT Auxiliary Structure ####
The manipulation of the parameters is managed by means of a structure variable named ``SAT`` which is created by the GUI and edited at the workspace or can also be created directly from workspace without the help of the GUI. This variable ``SAT`` has the following structure:

	> SAT
		SAT.NPIX								variable numeric
		SAT.name								array cell of string
		SAT.value								array cell of string or numeric
		SAT.attribute						array cell of strings

*	``.NPIX`` is a scalar structure member of ``SAT`` with the information of the number of footprints ``SAT`` is containing, for instance ``SAT.NPIX=369;`` indicates that the structure contains 369 footprints.
*	``.name`` is a cell structure member containing the variable names as strings, e.g. ``SAT.name{1}='SATELLITE_name'``, ``SAT.name{2}='Orbit_altitude_km'``, etc.
* ``.value`` is a cell structure member containing the values for the variable names ``.name``, for instance ``SAT.value{1}='GPM'``, ``SAT.value{2}=407``, an so on for the respectively examples in ``.name``.
* ``.attribute`` is a cell structure member with information about the variable in that specific cell, generally it contains a cell array with three inputs i.e. attribute type (global or var), variable type (number or string), and units ('m', or 'deg', etc.). As example ``SAT.name{2}`` has the attributes of ``SAT.attribute{2}={'global','number','km'};`` indicating that the ``Orbit_altitude_km`` variable has a global attribute with numeric value and unit in kilometer.

All information displayed in the [GUI](#edit-satellite-info-input) is then stored into the ``SAT`` structure variable which is used be the other scripts too.

### Read Satellite Info Input ###
 ``retrieve_SATinputdata.m`` This script encompass a function used to load the NetCDF satellite information file, it is used by the GUI ``edit_info_SAT4CMEM`` or can be run independently from workspace. The function runs without parameter or one or two parameters. When invoked without parameters then a Selection File GUI pops-up in order to browse a NetCDF file to open. When used with one parameter, that must be a string indicating the input file, and when used with two parameters then the first is a string with the input file and the second is a string with the absolute path to the input file.

 	USAGE:
 		> SAT = retrieve_SATinputdata;
		> SAT = retrieve_SATinputdata('input_netcdf_file.nc');
		> SAT = retrieve_SATinputdata('input_netcdf_file.nc','/path_to/file/');

It is mandatory to assign an output variable which is a structure of cell containing the data of the different variables of the satellite information file (``SAT`` see [example above](#sat-auxiliary-structure)).

### Create Satellite Info Input ###
 ``create_SATinput_nc.m`` This function is also called by ``edit_info_SAT4CMEM`` when the GUI introduced parameters want to be stored in a NetCDF file, but can also be used separately from the workspace.

When run without parameters (or by the GUI), the function first asks where to create the file and the name of it, once selected it creates a NetCDF file alike the ``SMOS_L1orbit_neckar.nc`` which can be used in CMEM by specifying the input parameter ``INPUTSATINFO`` (see section [Configuration Parameters](#cmem-configuration-parameters)).

The function also accept two parameters:

	> status = create_SATinput_nc(fname,SAT);

where ``fname`` is a string indicating the full path and the NetCDF file name, and the second argument is ``SAT`` which is the satellite structure containing all the information as indicated in section [Editing GUI](#sat-auxiliary-structure).

### Plotting CMEM Results ###

``plot_SATOPE_level4.m``: is a simple GUI script that allows to plot all the information included in a Lever 4 CMEM output, i.e. satellite simulator results for the input parameters assigned in the NetCDF satellite info at ``SATOPE_CLM4CMEM/forcing`` files. The Level 4 data is stored as a NetCDF file with the Brightness Temperatures (TB) at horizontal and vertical polarization for every pixel indicated with a longitude and latitude coordinates. This script plots TBs for every incidence angle presented in the NetCDF file as indicated by the variable ``THETA_INC``, and example is shown in Figure XX.

``plot_SATOPE_level1.m``: When the corresponding Level 1 data, related to the Level 4 data, is also available then that Level 1 data can also be plotted in a separated window. By clicking the button **Level 1** the script searches for a similar NetCDF files within the same directory with the flag ``lev1`` instead of ``lev4`` and plots them for all available incidence angle. Figure xx show the corresponding Level 1 data for the simulations in Figure xxx.

to do description of functions to plot and include Figures for the GNU Octave/ MATLAB scripts ...
