! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE RDCMEMGRIBAPIINFO 

! Purpose :
! -------
!     Get dimensions of GRIBAPI input data for cmem, and check their consistency
   
! Externals :
! ---------

! Internal datafields :
! -------------------
  
! Author :
!     Patricia de Rosnay, ECMWF, December 2008
!     Patricia de Rosnay  19-01-2012 add HTessel confi with variable LAI
!    P de Rosnay, ECMWF, May 2013, update netcdf and grib IO for multi-layer
!     
!------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM
USE YOMLUN   , ONLY : NULOUT
USE YOMCMEMPAR, ONLY : CIDVEG, CITVEG, CIATM, LGPRINT,CNAMEID, nlay_soil_ls
USE YOMCMEMGRIBAPI, ONLY : numberOfPoints,latitudeOfFirstgridPointInDegrees, &
      & longitudeOfFirstgridPointInDegrees,latitudeOfLastgridPointInDegrees,  &
      & longitudeOfLastgridPointInDegrees,CCVARAPI_NAME 
USE YOMCMEMFIELDS, ONLY : N, JJ, CLNAME 

IMPLICIT NONE



LOGICAL :: LD_firstcall
INTEGER(KIND=JPIM) :: i
CHARACTER(len=80) :: str_index

!
!------------------------------------------------------------------------------
!
! Read input files and get their dimensions   
!
!
!1- RD: Geopontential at surface
!-------------------------------
!
  CLNAME='Z.grib'
  CCVARAPI_NAME='Z'
  LD_firstcall = .True.
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  N = numberOfPoints
  LD_firstcall = .False.
!
!
!2- RD: snow depth
!-----------------
! 
  CLNAME='SD.grib'
  CCVARAPI_NAME='SD'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
! 
!3- RD: snow density
!-------------------
!   
  CLNAME='RSN.grib'
  CCVARAPI_NAME='RSN'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
!
!4- RD: Soil Temperature
!-----------------------
!
  DO i= 1, nlay_soil_ls
    write( str_index, * ) i
    str_index = adjustl(str_index)
    CLNAME='STL'//trim(str_index)//'.grib'
    CCVARAPI_NAME='STL'//trim(str_index)
    CALL io_cmemgribapi(CLNAME,LD_firstcall)
    IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
  END DO
!  
  CLNAME='TSKIN.grib'
  CCVARAPI_NAME='TSKIN'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
!
!5- RD: Soil moisture
!--------------------
!
  DO i= 1, nlay_soil_ls
    write( str_index, * ) i
    str_index = adjustl(str_index)
    CLNAME='SWVL'//trim(str_index)//'.grib'
    CCVARAPI_NAME='SWVL'//trim(str_index)
    CALL io_cmemgribapi(CLNAME,LD_firstcall)
    IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
  END DO

!
!      
!6- Soil Texture
!---------------
!
  CLNAME='sand.grib'
  CCVARAPI_NAME='SAND'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
  CLNAME='clay.grib'
  CCVARAPI_NAME='CLAY'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
!
!7- RD: Vegetation 
!-----------------
!   
SELECT CASE (CIDVEG)
    
 CASE ( 'Ecoclimap' )  
!
! low vegetation fraction
  CLNAME='ECOCVL.grib'
  CCVARAPI_NAME='CVL'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
! high vegetation fraction
  CLNAME='ECOCVH.grib'
  CCVARAPI_NAME='CVH'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
! low vegetation types 
  CLNAME='ECOTVL.grib'
  CCVARAPI_NAME='TVL'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
! high vegetation types 
  CLNAME='ECOTVH.grib'
  CCVARAPI_NAME='TVH'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
! water fraction
  CLNAME='ECOWAT.grib'
  CCVARAPI_NAME='WATER'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
! LAI of low vegetation 
  CLNAME='ECOLAIL.grib'
  CCVARAPI_NAME='LAI'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)

 CASE ('Tessel','HTessel' )

! low vegetation fraction
  CLNAME='CVL.grib'
  CCVARAPI_NAME='CVL'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
! high vegetation fraction
  CLNAME='CVH.grib'
  CCVARAPI_NAME='CVH'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
! low vegetation type
  CLNAME='TVL.grib'
  CCVARAPI_NAME='TVL'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
! high vegetation type
  CLNAME='TVH.grib'
  CCVARAPI_NAME='TVH'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
!
! water fraction 
  CLNAME='LSM.grib' 
  CCVARAPI_NAME='LSM'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)

  IF ( CIDVEG == 'HTessel') THEN             
    CLNAME='LAIL.grib'
    CCVARAPI_NAME='LAI'
    CALL io_cmemgribapi(CLNAME,LD_firstcall)
    IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
  ENDIF

!
END SELECT
!
!
!8- RD:  Vegetation and air temperature
!---------------------------------------
!
SELECT CASE (CITVEG)
  CASE ( 'Tair' )
  ! RD: T air 2m
  CLNAME='2T.grib'
  CCVARAPI_NAME='TAIR'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
END SELECT
!
SELECT CASE (CIATM)
  CASE ( 'Pellarin', 'Ulaby' )
  ! RD: T air 2m
  CLNAME='2T.grib'
  CCVARAPI_NAME='TAIR'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  IF(N .ne. numberOfPoints) CALL ABOR1('Non consistent file size'//CLNAME)
END SELECT
!
!  

WRITE(NULOUT,*) 'Forcing_info: input files are consistent.'
WRITE(NULOUT,*) 'Number of points: ',numberOfPoints
WRITE(NULOUT,*) 'Longitude range: ',longitudeOfFirstGridPointInDegrees,longitudeOfLastGridPointInDegrees
WRITE(NULOUT,*) 'Latitude range: ',latitudeOfFirstGridPointInDegrees,latitudeOfLastGridPointInDegrees

         
END SUBROUTINE RDCMEMGRIBAPIINFO
