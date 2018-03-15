! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE RDCMEMGRIBAPI

! Purpose :
! -------
!     Read input data for cmem from GRIBAPI files
   
! Externals :
! ---------

! Internal datafields :
! -------------------
  
! Author :
!     Patricia de Rosnay, ECMWF, December 2008
!     Patricia de Rosnay, ECMWF, November 2009 SSS 
!     Patricia de Rosnay, ECMWF, 19-01-2012 add HTessel confi with variable LAI
!    P de Rosnay, ECMWF, May 2013, update netcdf and grib IO for multi-layer
!
!     
!------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPRM
USE YOMLUN   , ONLY : NULOUT
USE YOMCMEMPAR, ONLY : CIDVEG, CITVEG, CIATM, LGPRINT,CNAMEID, PPG &
                    & , LOMASK_OCEAN,LOMASK_AUTO &
                     & , nlay_soil_ls
USE YOMCMEMGRIBAPI, ONLY : fapifield, numberOfValues,CCVARAPI_NAME
USE YOMCMEMFIELDS, ONLY : N, JJ, CLNAME, fTVL, fTVH,  fs_laiL, ftfrac, fs_laiL &
                    &  , fsnowd, frsnow, fwater, ftl_lsm, ftair,ftveg &
                    & , ftskin,fwc_lsm,fsand,fclay,fwater, mask 
USE YOMCMEMATM, ONLY : fZ
USE YOMCMEMSOIL, ONLY : sal_sea

IMPLICIT NONE

INTEGER(KIND=JPIM) :: i
REAL(KIND=JPRM) :: rvcov(0:20)
REAL(KIND=JPRM) :: rvlai(0:20)
REAL(KIND=JPRM) :: b(0:7)
REAL(KIND=JPRM) :: b1(0:7)
REAL(KIND=JPRM) :: b2(0:7)
REAL(KIND=JPRM) :: b3(0:7)
REAL(KIND=JPRM) :: VWC(0:7)
REAL(KIND=JPRM) :: ZNrh(0:7)
REAL(KIND=JPRM) :: ZNrv(0:7)
REAL(KIND=JPRM) :: Ztth(0:7)
REAL(KIND=JPRM) :: Zttv(0:7)
REAL(KIND=JPRM) :: Zhr(0:7)
REAL(KIND=JPRM) :: Zw_eff(0:7,2)
INTEGER(KIND=JPIM) :: RVTV(0:20)

REAL(KIND=JPRM),ALLOCATABLE :: cvegl(:), cvegh(:)
REAL(KIND=JPRM),ALLOCATABLE :: fbare(:)  
REAL(KIND=JPRM),ALLOCATABLE :: fvegl(:)  
REAL(KIND=JPRM),ALLOCATABLE :: fvegh(:)    
REAL(KIND=JPRM),ALLOCATABLE :: sncov(:)

INTEGER(KIND=JPIM) :: JJ2D,JTIME,N2D
LOGICAL :: LD_firstcall
CHARACTER(len=80) :: str_index

!
!------------------------------------------------------------------------------
!
! Read input files 
!
!
  LD_firstcall = .False.
!
!Define default mask 
!
mask(:) = 1_JPIM
!
!
!1- RD: Geopontential at surface (in km)
!---------------------------------------
!
  CLNAME='Z.grib'
  CCVARAPI_NAME='Z'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  fZ(:) = fapifield(:,1)
!
  DO JJ = 1, N
    fZ(JJ)=fZ(JJ)/1000._JPRM/PPG
  ENDDO
!
! Update the mask to eliminate unrealistic points
 SELECT CASE (LOMASK_AUTO)
     CASE (.True.)
     WHERE (fZ > 9.0_JPRM .or. fZ < -1.0_JPRM ) mask = 0_JPIM  ! unexpected elevation on the Earth
 END SELECT
!
!
!2- RD: snow depth (m eqv water)
!-------------------------------
! 
  CLNAME='SD.grib'
  CCVARAPI_NAME='SD'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  fsnowd(:) = fapifield(:,1)
  ! Minimum thresold of snow to be considered
  ! To avoid numerical escillation of HUT model in case of thin snow layer:
  fsnowd(:) = MAX(1.e-2_JPRM,fsnowd(:))
  WHERE(fsnowd == 1.e-2_JPRM) fsnowd(:) = 0.0
!
! 
!3- RD: snow density (kg/m3)
!--------------------------
!   
  CLNAME='RSN.grib'
  CCVARAPI_NAME='RSN'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  frsnow(:) = fapifield(:,1)
!
!
!4- RD: Soil Temperature (K)
!---------------------------
!
  DO i= 1, nlay_soil_ls
    write( str_index, * ) i
    str_index = adjustl(str_index)
    CLNAME='STL'//trim(str_index)//'.grib'
    CCVARAPI_NAME='STL'//trim(str_index)
    CALL io_cmemgribapi(CLNAME,LD_firstcall)
    ftl_lsm(:,i) = fapifield(:,1)
  END DO

!
  CLNAME='TSKIN.grib'
  CCVARAPI_NAME='TSKIN'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  ftskin(:) = fapifield(:,1)
!
! Update the mask to eliminate unrealistic points
 SELECT CASE (LOMASK_AUTO)
  CASE (.True.)
  do i = 1,nlay_soil_ls
    WHERE (ftl_lsm(:,i) <100.0_JPRM .or.  ftskin < 100.0_JPRM) mask(:) = 0_JPIM  ! unrealistic temperature
  enddo
  END SELECT
!
!
!5- RD: Soil moisture (m3/m3)
!----------------------------
!
  DO i= 1, nlay_soil_ls
    write( str_index, * ) i
    str_index = adjustl(str_index)
    CLNAME='SWVL'//trim(str_index)//'.grib'
    CCVARAPI_NAME='SWVL'//trim(str_index)
    CALL io_cmemgribapi(CLNAME,LD_firstcall)
    fwc_lsm(:,i) = fapifield(:,1)
  END DO
!
! Update the mask to eliminate unrealistic points
 SELECT CASE (LOMASK_AUTO)
   CASE (.True.)
   do i = 1, nlay_soil_ls
     WHERE (fwc_lsm(:,i) < 0.0_JPRM ) mask = 0_JPIM  ! Unrealistic SSM 
   enddo
 END SELECT

   ! in any case ensure sm is positive
   fwc_lsm = MAX(0.0_JPRM, fwc_lsm)
!
!      
!6- RD: Soil Texture (%)
!-----------------------
!
  CLNAME='sand.grib'
  CCVARAPI_NAME='SAND'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
    fsand(:) = fapifield(:,1)
!
  CLNAME='clay.grib'
  CCVARAPI_NAME='CLAY'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
    fclay(:) = fapifield(:,1)
!
!
!7- RD: Vegetation 
!-----------------
!   
! ALLOCATE TEMPORARY VARIABLES
  ALLOCATE (fbare(N))
  ALLOCATE (fvegl(N))
  ALLOCATE (fvegh(N))
  ALLOCATE (cvegl(N))
  ALLOCATE (cvegh(N))
  ALLOCATE (sncov(N))

  cvegh(:) = 1.0_JPRM
  cvegl(:) = 1.0_JPRM

SELECT CASE (CIDVEG)
    
 CASE ( 'Ecoclimap' )  
!
! low vegetation fraction (-)
  CLNAME='ECOCVL.grib'
  CCVARAPI_NAME='CVL'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
    fvegl(:) = fapifield(:,1)
!
! high vegetation fraction (-)
  CLNAME='ECOCVH.grib'
  CCVARAPI_NAME='CVH'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
    fvegh(:) = fapifield(:,1)
!
! low vegetation types 
  CLNAME='ECOTVL.grib'
  CCVARAPI_NAME='TVL'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
    fTVL(:) = fapifield(:,1)
!
! high vegetation types 
  CLNAME='ECOTVH.grib'
  CCVARAPI_NAME='TVH'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
    fTVH(:) = fapifield(:,1)
!
! water fraction (-)
  CLNAME='ECOWAT.grib'
  CCVARAPI_NAME='WATER'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
    fwater(:) = fapifield(:,1)

!
! LAI of low vegetation 
  CLNAME='ECOLAIL.grib'
  CCVARAPI_NAME='LAI'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  fs_laiL(:) = fapifield(:,1)

!Clean up vegetation tables
  WHERE ( fTVL > 7_JPIM .and. fwater(:) <  0.5_JPRM ) fTVL = 4_JPIM  ! Defaults low vegetation is grass
  WHERE ( fTVL > 7_JPIM .and. fwater(:) >= 0.5_JPRM ) fTVL = 0_JPIM  ! No vegetation on  water pixel
  WHERE ( fTVH > 7_JPIM .and. fwater(:) <  0.5_JPRM ) fTVH = 1_JPIM  ! Defaults low vegetation is grass
  WHERE ( fTVH > 7_JPIM .and. fwater(:) >= 0.5_JPRM ) fTVH = 0_JPIM  ! No vegetation on  water pixel


 CASE ('Tessel','HTessel' )

! low vegetation fraction (-)
  CLNAME='CVL.grib'
  CCVARAPI_NAME='CVL'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
    fvegl(:) = fapifield(:,1)
!
! high vegetation fraction (-)
  CLNAME='CVH.grib'
  CCVARAPI_NAME='CVH'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
    fvegh(:) = fapifield(:,1)
!
! low vegetation type
  CLNAME='TVL.grib'
  CCVARAPI_NAME='TVL'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
    fTVL(:) = fapifield(:,1)

!
! high vegetation type
  CLNAME='TVH.grib'
  CCVARAPI_NAME='TVH'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
    fTVH(:) = fapifield(:,1)
!
! water fraction (from LSM in TESSEL)  (-)
  CLNAME='LSM.grib' 
  CCVARAPI_NAME='LSM'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
    fwater(:) =  1.0_JPIM - fapifield(:,1)
!
    ! LAI of low vegetation 
    !
    CALL VEGTABLE (RVLAI,RVCOV,RVTV,b,VWC,b1,b2,b3,ZNrh,ZNrv,Ztth,Zttv,Zhr,Zw_eff)
    IF ( CIDVEG == 'Tessel') THEN             
      fs_laiL(:) = rvlai(fTVL(:))
    ELSE 
    ! LAI of low vegetation 
      CLNAME='LAIL.grib'
      CCVARAPI_NAME='LAI'
     CALL io_cmemgribapi(CLNAME,LD_firstcall)
     fs_laiL(:) = fapifield(:,1)
    ENDIF
END SELECT
!
!7.1- Calculate tile fractions
!   
    ! TESSEL vegetation map needs to account for bare soil
    
    IF ( CIDVEG == 'Tessel' .OR. CIDVEG == 'HTessel') THEN             
      cvegl(:) = rvcov(fTVL(:))  
      cvegh(:) = rvcov(fTVH(:))
      ! Convert model vegetation types to CMEM vegetation requirement
      fTVL(:) = RVTV(fTVL(:))
      fTVH(:) = RVTV(fTVH(:))
    ENDIF

  
    fvegl(:) =  fvegl(:) * cvegl(:)
    fvegh(:) =  fvegh(:) * cvegh(:)
    fbare(:) = 1. - fvegl(:) - fvegh(:)

    ! Snow cover on land tiles

    WHERE (  fsnowd(:) > 0.0_JPRM) sncov(:) = 1.0_JPRM
    WHERE (  fsnowd(:) <= 0.0_JPRM) sncov(:) = 0.0_JPRM


    ftfrac(:,1) = fbare(:) * (1. - fwater(:)) * (1.-sncov(:)) ! bare soil
    ftfrac(:,2) = fbare(:) * (1. - fwater(:)) * sncov(:)      ! bare soil with snow
    ftfrac(:,3) = fvegl(:) * (1. - fwater(:)) * (1.-sncov(:)) ! low veg
    ftfrac(:,4) = fvegl(:) * (1. - fwater(:)) * sncov(:)      ! low veg with snow
    ftfrac(:,5) = fvegh(:) * (1. - fwater(:)) * (1.-sncov(:)) ! high veg
    ftfrac(:,6) = fvegh(:) * (1. - fwater(:)) * sncov(:)      ! high veg with snow
    ftfrac(:,7) = fwater(:)                                ! water

    SELECT CASE (LOMASK_OCEAN)
       CASE (.True.)
        WHERE ( fwater(:) >= 0.5 )  mask(:) = 0.0_JPRM
    END SELECT
  




! DEALLOCATE TEMPORARY VARIABES
  DEALLOCATE (fbare)
  DEALLOCATE (fvegl)
  DEALLOCATE (fvegh)
  DEALLOCATE (cvegl)
  DEALLOCATE (cvegh)
  DEALLOCATE (sncov)
!
!
!8- RD:  Vegetation and air temperatures
!---------------------------------------
!
SELECT CASE (CITVEG)
  CASE ( 'Tsurf' )
  ftveg(:) = ftl_lsm(:,1)
  CASE ( 'Tair' )
  ! RD: T air 2m (K)
  CLNAME='2T.grib'
  CCVARAPI_NAME='TAIR'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  ftveg(:) = fapifield(:,1)
!
! Update the mask to eliminate unrealistic points
SELECT CASE (LOMASK_AUTO)
    CASE (.True.)
    WHERE (ftveg < 100.0_JPRM ) mask = 0_JPIM  ! Unrealistic Temperature 
END SELECT
!
END SELECT
!
SELECT CASE (CIATM)
  CASE ( 'Pellarin', 'Ulaby' )
  ! RD: T air 2m (K)
  CLNAME='2T.grib'
  CCVARAPI_NAME='TAIR'
  CALL io_cmemgribapi(CLNAME,LD_firstcall)
  ftair(:) = fapifield(:,1)
!
! Update the mask to eliminate unrealistic points
  SELECT CASE (LOMASK_AUTO)
     CASE (.True.)
     WHERE (ftair < 100.0_JPRM ) mask = 0_JPIM  ! Unrealistic Temperature 
  END SELECT
!
ENDSELECT
!
!  
!-------------------------------------------------------------------
! 9.  Read SSS
!-------------------------------------------------------------------

sal_sea(:) = 32.5_JPRM

         
END SUBROUTINE RDCMEMGRIBAPI
