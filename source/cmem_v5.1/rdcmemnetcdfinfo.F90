! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE RDCMEMNETCDFINFO 

! Purpose :
! -------
!     Get dimensions of netcdf input data for cmem, and check their consistency
   
! Externals :
! ---------

! Internal datafields :
! -------------------
  
! Author :
!     Patricia de Rosnay, ECMWF, January 2008
!                         19-01-2012 add HTessel confi with variable LAI
!    P de Rosnay, ECMWF, May 2013, update netcdf and grib IO for multi-layer
!     
!------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM
USE YOMLUN   , ONLY : NULOUT
USE YOMCMEMPAR, ONLY : CIDVEG, CITVEG, CIATM, LGPRINT,CNAMEID, nlay_soil_ls
USE YOMCMEMNETCDF, ONLY : NTIMES, NDIMS, NLATS, NLONS, NLVLS, CCVARNC_NAME &
                    &   , VAR_IN,xlats,xlons,xtimes,xlvls
USE YOMCMEMFIELDS, ONLY : N, JJ, CLNAME 

IMPLICIT NONE



INTEGER(KIND=JPIM) :: JJ2D,JTIME,N2D,i
LOGICAL :: LD_firstcall, LD_infocoord
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
  CLNAME='../forcing/Z.nc'
  CCVARNC_NAME='Z'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  N2D = NLONS*NLATS
!
!
!2- RD: snow depth
!-----------------
! 
  CLNAME='../forcing/SD.nc'
  CCVARNC_NAME='SD'
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size')
  N = NLONS*NLATS*NTIMES
!
! 
!3- RD: snow density
!-------------------
!   
  CLNAME='../forcing/RSN.nc'
  CCVARNC_NAME='RSN'
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
  IF(N .ne. NLONS*NLATS*NTIMES) CALL ABOR1('Non consistent file size N'//CLNAME)
!
!
!4- RD: Soil Temperature
!-----------------------
!
  DO i= 1, nlay_soil_ls
    write( str_index, * ) i
    str_index = adjustl(str_index)
    CLNAME='../forcing/STL'//trim(str_index)//'.nc'
    CCVARNC_NAME='STL'//trim(str_index)
    CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
    IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
    IF(N .ne. NLONS*NLATS*NTIMES) CALL ABOR1('Non consistent file size N'//CLNAME)
  END DO

!  
  CLNAME='../forcing/TSKIN.nc'
  CCVARNC_NAME='TSKIN'
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
  IF(N .ne. NLONS*NLATS*NTIMES) CALL ABOR1('Non consistent file size N'//CLNAME)
!
!
!5- RD: Soil moisture
!--------------------
!
  DO i= 1, nlay_soil_ls
    write( str_index, * ) i
    str_index = adjustl(str_index)
    CLNAME='../forcing/SWVL'//trim(str_index)//'.nc'
    CCVARNC_NAME='SWVL'//trim(str_index)
    LD_firstcall = .True.
    LD_infocoord = .True.
    CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
    IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
    IF(N .ne. NLONS*NLATS*NTIMES) CALL ABOR1('Non consistent file size N'//CLNAME)
  END DO

!
!      
!6- Soil Texture
!---------------
!
  CLNAME='../forcing/sand.nc'
  CCVARNC_NAME='SAND'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
!
  CLNAME='../forcing/clay.nc'
  CCVARNC_NAME='CLAY'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
!
!
!7- RD: Vegetation 
!-----------------
!   
SELECT CASE (CIDVEG)
    
 CASE ( 'Ecoclimap' )  
!
! low vegetation fraction
  CLNAME='../forcing/ECOCVL.nc'
  CCVARNC_NAME='CVL'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
!
! high vegetation fraction
  CLNAME='../forcing/ECOCVH.nc'
  CCVARNC_NAME='CVH'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
!
! low vegetation types 
  CLNAME='../forcing/ECOTVL.nc'
  CCVARNC_NAME='TVL'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
!
! high vegetation types 
  CLNAME='../forcing/ECOTVH.nc'
  CCVARNC_NAME='TVH'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
!
! water fraction
  CLNAME='../forcing/ECOWAT.nc'
  CCVARNC_NAME='WATER'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
!
! LAI of low vegetation 
  CLNAME='../forcing/ECOLAIL.nc'
  CCVARNC_NAME='LAI'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)

 CASE ('Tessel','HTessel' )

! low vegetation fraction
  CLNAME='../forcing/CVL.nc'
  CCVARNC_NAME='CVL'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
!
! high vegetation fraction
  CLNAME='../forcing/CVH.nc'
  CCVARNC_NAME='CVH'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
!
! low vegetation type
  CLNAME='../forcing/TVL.nc'
  CCVARNC_NAME='TVL'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
!
! high vegetation type
  CLNAME='../forcing/TVH.nc'
  CCVARNC_NAME='TVH'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
!
! water fraction 
  CLNAME='../forcing/LSM.nc' 
  CCVARNC_NAME='LSM'
  LD_firstcall = .False.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
!
  IF ( CIDVEG == 'HTessel') THEN             
    CLNAME='../forcing/LAIL.nc'
    CCVARNC_NAME='LAI'
    LD_firstcall = .True.
    LD_infocoord = .False.
    CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
    IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
  ENDIF
END SELECT
!
!
!8- RD:  Vegetation and air temperature
!---------------------------------------
!
SELECT CASE (CITVEG)
  CASE ( 'Tair' )
  ! RD: T air 2m
  CLNAME='../forcing/2T.nc'
  CCVARNC_NAME='TAIR'
  LD_firstcall = .False.
  LD_infocoord = .True.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
  IF(N .ne. NLONS*NLATS*NTIMES) CALL ABOR1('Non consistent file size N'//CLNAME)
END SELECT
!
SELECT CASE (CIATM)
  CASE ( 'Pellarin', 'Ulaby' )
  ! RD: T air 2m
  CLNAME='../forcing/2T.nc'
  CCVARNC_NAME='TAIR'
  LD_firstcall = .True.
  LD_infocoord = .False.
  CALL io_cmemnetcdf(CLNAME,LD_firstcall,LD_infocoord)
  IF(N2D .ne. NLONS*NLATS) CALL ABOR1('Non consistent file size'//CLNAME)
  IF(N .ne. NLONS*NLATS*NTIMES) CALL ABOR1('Non consistent file size N'//CLNAME)
END SELECT
!
!  

WRITE(NULOUT,*) 'Forcing_info: input files are consistent.'
WRITE(NULOUT,*) 'Lat x Lon Grid dimension: ',N2D 
WRITE(NULOUT,*) 'Lat x Lon x Time Grid dimension: ',N 
WRITE(NULOUT,*) 'Longitude range: ',xlons(1),xlons(NLONS) 
WRITE(NULOUT,*) 'Latitude range: ',xlats(1),xlats(NLATS) 
WRITE(NULOUT,*) 'Time range: ',1,NTIMES 

         
END SUBROUTINE RDCMEMNETCDFINFO
