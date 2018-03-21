! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
PROGRAM CMEM_MAIN
!
!
!
!       PURPOSE.
!       -------
!     CMEM Main program
!     CMEM: Copyright © ECMWF
!     Compute Low frequency passive microwave emission of the surface,
!     from 1.4 to 20 GHz.
!
!       INTERFACE :
!       ---------
!
!       EXTERNALS :
!       ---------
!         RD_*    : Read input data
!         CMEM_*  : module of emission model
!         IO_* : read grib and netcdf data
!         WR*  : write output data
!
!       REFERENCES
!      ------------
!      Drusch, M., T. Holmes, P. de Rosnay and G. Balsamo, 2009:
!      Comparing ERA-40 based L-band brightness temperatures with Skylab observations:
!      A calibration / validation study using the
!      Community Microwave Emission Model, J. Hydromet., Vol 10, DOI: 10.1175/2008JHM964.1, 2009
!
!      de Rosnay P., M. Drusch, A. Boone, G. Balsamo,
!      B. Decharme, P. Harris, Y. Kerr, T. Pellarin, J.
!      Polcher and J.-P. Wigneron,
!      "The AMMA Land Surface Model Intercomparison
!      Experiment coupled to the Community Microwave Emission Model: ALMIP-MEM",
!      J. Geophys. Res., Vol 114, doi:10.1029/2008JD010724, 2009
!
!      AUTHORS.
!      --------
!     Thomas Holmes   ECMWF     v0 6/09/06
!     Patricia de Rosnay ECMWF: v1.1 Dec  2007 (Grib and ASCII options)
!                               v1.2 February 2008 (Netcdf option added)
!                               v1.3 April 2008
!                               v1.4 November 2008
!                                    December 2008 (Grib API option added)
!                               v2.0 January 2009  (GRIB API version released)
!                               v2.1 March 2009  (gfortran and IFC compatibility)
!                               v3.0 November 2009  (IFS interface structure, SSS(N))
!                               v4.0 February 2012  HTessel option, bug fixes in TEFF, Liebe and VEGTABLE
!                               v4.1 May 2012  Improved Wilheit multilayer (flexible LSM and MW layers)
!                               v5.1 February 2014 fixed bugs, code cleaning and Apache license
!
!------------------------------------------------------------------------------
! ===================================================================
! CODE HAS BEEN MODIFIED FROM THE ORIGINAL CMEM_MAIN.F90 SOURCE CODE
! this subroutine allows to call CMEM with the SATELLITE OPERATOR from
! any other source code needed.
!
! Modified by:
! (c) 2017 P. Saavedra Garfias (pablosaa@uni-bonn.de) UNI BONN
! see: LICENCE
! -------------------------------------------------------------------

!
USE PARKIND1, ONLY : JPIM, JPRM
USE YOMLUN, ONLY : NULOUT, NULTMP
USE YOMCMEMPAR
USE YOMCMEMFIELDS
USE YOMCMEMSOIL
USE YOMCMEMVEG, ONLY : wc_veg, tb_veg, w_eff, bj, tauN, &
                      & tb_veg, tau_veg, t_veg, a_geo, a_geoL, a_geoH,tth,ttv
USE YOMCMEMATM, ONLY : tau_atm, tb_au, tb_ad, tb_toa, tb_tov,fZ,tair
use rdclm_wrcmem ! PSG: module with subroutines to read/write CLM NetCDF
use clm4cmem ! PSG: module with toolbox for CLM integration
use yomcmemnetcdf, only: NTIMES, NLATS, NLONS, NINC  ! PSG: compativility
use SatOperator ! PSG: module with toolbox for Satellite Operator
!
IMPLICIT NONE
!
INTEGER(KIND=JPIM) :: JJPOL
REAL(KIND = JPRM) :: RSN(2), ESN(2)
REAL(KIND = JPRM) :: tfrac(JPCMEMTILE)
!
!
! 1.0 Model SETUP
! ------------------------------------
integer(kind=JPIM) :: Nvars(3)        ! PSG:
type(CLM_DATA) :: CLMVARS             ! PSG:
type(SATELLITE) :: SAT   ! PSG:
! Assigning input parameters to the CMEM global parameters
call getarg(1,CLNAME)
call getarg(2,INPUTNAMLST)
print*, "CLM file: ",CLNAME
!
INCLUDE "core_cmemcode.F90"
!

!
deallocate(SAT%TBSAT_HV,SAT%lon_foprt,SAT%lat_foprt,SAT%incl_foprt,SAT%theta)
END PROGRAM CMEM_MAIN

! end of main CMEM Stand-alone program.
