! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE RDCMEMASCII 

! Purpose :
! -------
!     Read input data for cmem from grib files
   
! Externals :
! ---------

! Internal datafields :
! -------------------
  
!     Patricia de Rosnay  01-10-2007
!     Patricia de Rosnay  02-11-2009 SSS
!     Patricia de Rosnay  February 2011 fix in multi-layer atmospheric model
!     Patricia de Rosnay, 19-01-2012 add HTessel confi with variable LAI
!------------------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRM
USE YOMLUN, ONLY : NULOUT, NULTMP
USE YOMCMEMPAR, ONLY : CIDVEG, CITVEG, CIATM,  LGPRINT,CNAMEID, LOFIELDEXP, PPG, &
                    & LOMASK_OCEAN, LOMASK_AUTO,nlay_soil_ls
USE YOMCMEMFIELDS, ONLY : N, JJ, doy,ndate,ntime,nbpt &
                    &  ,fTVL, fTVH,  fs_laiL, ftfrac, fsnowd, frsnow,ftl_lsm, ftair,ftveg &
                    &  ,ftskin,fwc_lsm,fsand,fclay,fwater,mask 
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

REAL(KIND=JPRM) ::  zdoy,zsand,zclay, zfZ,zfTVL,zfTVH,zfvegl,zfvegh,zfwater,zlaiH

!
!------------------------------------------------------------------------------
!
!
!Define default mask
!
mask(:) = 1_JPIM
!



!----------------------------------------------------------
! 1.  Read main input file 
!----------------------------------------------------------

OPEN(NULTMP,FILE='forcing_cmem_main.asc',status='old')



DO JJ = 1,N
READ (NULTMP,*) nbpt(JJ), ndate(JJ),ntime(JJ),(fwc_lsm(JJ,i),i=1,nlay_soil_ls),ftskin(JJ)&
              & ,(ftl_lsm(JJ,i),i=1,nlay_soil_ls),fsnowd(JJ),frsnow(JJ),doy(JJ)
ENDDO 
CLOSE(NULTMP)


fsnowd(:) = MAX(1.e-2_JPRM,fsnowd(:))
WHERE(fsnowd == 1.e-2_JPRM) fsnowd(:) = 0.0

! Update the mask to eliminate unrealistic points
   SELECT CASE (LOMASK_AUTO)
   CASE (.True.)
     do i = 1,nlay_soil_ls
      WHERE (ftl_lsm(:,i) <100.0_JPRM .or. ftskin(:) < 100.0_JPRM) mask(:) = 0_JPIM  ! unrealistic temperature
      WHERE ( fwc_lsm(:,i) < 0.0_JPRM ) mask(:) = 0_JPIM  ! Unrealistic SM
     enddo
   END SELECT

     ! in any case ensure sm is positive
     fwc_lsm(:,:) = MAX(0.0_JPRM, fwc_lsm(:,:))


!----------------------------------------------------------
! 2.  Read soil and vegetation input files
!----------------------------------------------------------

! 2.1 Read data texture and vegetation types
!---------------------------------------------

  ALLOCATE (fbare(N))
  ALLOCATE (fvegl(N))
  ALLOCATE (fvegh(N))
  ALLOCATE (cvegl(N))
  ALLOCATE (cvegh(N))
  ALLOCATE (sncov(N))

  ! default vegetation cover is 100 % on vegetated tiles
  cvegh(:) = 1.0_JPRM
  cvegl(:) = 1.0_JPRM



SELECT CASE (LOFIELDEXP)
       
   CASE (.TRUE.)  ! Field experiment

       ! constant parameters to be read (%, %, m2/s2)

       OPEN(NULTMP,FILE='forcing_cmem_soil_constant.asc',status='old')
       READ (NULTMP,*) zsand,zclay,zfZ
       CLOSE(NULTMP)
       ! Convert geopotential at surface from (m2/s2) to height (km):    
       zfZ = zfZ / 1000._JPRM / PPG

       ! Update the mask to eliminate unrealistic points
       SELECT CASE (LOMASK_AUTO)
         CASE (.True.)
         WHERE (fZ > 9.0_JPRM .or. fZ < -1.0_JPRM ) mask = 0_JPIM  ! unexpected elevation on the Earth
       END SELECT

       SELECT CASE (CIDVEG)

          CASE ('Tessel','HTessel')

           ! fraction in (-)
           OPEN(NULTMP,FILE='forcing_cmem_veg-tessel_constant.asc',status='old')
           READ (NULTMP,*) zfTVL,zfTVH,zfvegl,zfvegh,zfwater
           CLOSE(NULTMP)

          CASE ('Ecoclimap')

           ! fraction in (-)
           OPEN(NULTMP,FILE='forcing_cmem_veg-ecoclimap_constant.asc',status='old')
           WRITE(NULOUT,*) 'Read forcing_cmem_veg'
           READ (NULTMP,*) zfTVL,zfTVH,zfvegl,zfvegh,zfwater
           CLOSE(NULTMP)

       END SELECT

       ! Apply constant values to 1D file
       fsand = zsand
       fclay = zclay
       fZ = zfZ
       fTVL = zfTVL
       fTVH = zfTVH
       fvegl = zfvegl
       fvegh = zfvegh
       fwater = zfwater

  CASE (.FALSE.)

       ! Surface conditions vary
        OPEN(NULTMP,FILE='forcing_cmem_soil.asc',status='old')
        DO JJ = 1,N
         READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),fsand(JJ),fclay(JJ),fZ(JJ),doy(JJ)
        ENDDO
       ! Convert geopotential at surface from (m2/s2) to height (km):    
       fZ = fZ / 1000._JPRM / PPG
        CLOSE(NULTMP)

       SELECT CASE (CIDVEG)

         CASE ('Tessel','HTessel')

           OPEN(NULTMP,FILE='forcing_cmem_veg-tessel.asc',status='old')
           WRITE(NULOUT,*) 'Read forcing_cmem_veg-tessel'
           DO JJ = 1,N
             READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),fTVL(JJ),fTVH(JJ),fvegl(JJ),fvegh(JJ),fwater(JJ),doy(JJ)
           ENDDO
           CLOSE(NULTMP)

         CASE ('Ecoclimap')

           OPEN(NULTMP,FILE='forcing_cmem_veg-ecoclimap.asc',status='old')
           WRITE(NULOUT,*) 'Read forcing_cmem_veg-ecoclimap'

           DO JJ = 1,N
             READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),fTVL(JJ),fTVH(JJ),fvegl(JJ),fvegh(JJ),fwater(JJ),doy(JJ)
           ENDDO
           CLOSE(NULTMP)
         
           WHERE ( fTVL > 7_JPIM .and. fwater(:) <  0.5_JPRM ) fTVL = 4_JPIM  ! Defaults low vegetation is grass
           WHERE ( fTVL > 7_JPIM .and. fwater(:) >= 0.5_JPRM ) fTVL = 0_JPIM  ! No vegetation on  water pixel
           WHERE ( fTVH > 7_JPIM .and. fwater(:) <  0.5_JPRM ) fTVH = 1_JPIM  ! Defaults low vegetation is grass
           WHERE ( fTVH > 7_JPIM .and. fwater(:) >= 0.5_JPRM ) fTVH = 0_JPIM  ! No vegetation on  water pixel

       END SELECT


  

END SELECT  ! LOFIELDEXP


! 2.2  Read lai 
!--------------


SELECT CASE (CIDVEG)

   CASE ('Ecoclimap')

    OPEN(NULTMP,FILE='forcing_cmem_lai.asc',status='old')


    DO JJ = 1,N
    READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),fs_laiL(JJ),zlaiH,doy(JJ)
    ENDDO
    CLOSE(NULTMP)

    CASE ('Tessel')

    CALL VEGTABLE (RVLAI,RVCOV,RVTV,b,VWC,b1,b2,b3,ZNrh,ZNrv,Ztth,Zttv,Zhr,Zw_eff)
    fs_laiL(:) = rvlai(fTVL(:))
    ! TESSEL vegetation map needs to account for bare soil
    ! Once LAI is computed, convert to ECOCLIMAP
    cvegl(:) = rvcov(fTVL(:)) 
    cvegh(:) = rvcov(fTVH(:))  
    ! Convert model vegetation types to CMEM vegetation requirement
    fTVL(:) = RVTV(fTVL(:))
    fTVH(:) = RVTV(fTVH(:))

    CASE ('HTessel')
    ! read LAI
    OPEN(NULTMP,FILE='forcing_cmem_lai.asc',status='old')
    DO JJ = 1,N
    READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),fs_laiL(JJ),zlaiH,doy(JJ)
    ENDDO
    CLOSE(NULTMP)
    ! HTESSEL vegetation map needs to account for bare soil
    ! Once LAI is computed, convert to ECOCLIMAP
    CALL VEGTABLE (RVLAI,RVCOV,RVTV,b,VWC,b1,b2,b3,ZNrh,ZNrv,Ztth,Zttv,Zhr,Zw_eff)
    cvegl(:) = rvcov(fTVL(:)) 
    cvegh(:) = rvcov(fTVH(:))  
    ! Convert model vegetation types to CMEM vegetation requirement
    fTVL(:) = RVTV(fTVL(:))
    fTVH(:) = RVTV(fTVH(:))



END SELECT


! 2.3  Calculate tile fractions
!------------------------------

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



  DEALLOCATE (fbare)
  DEALLOCATE (fvegl)
  DEALLOCATE (fvegh)
  DEALLOCATE (cvegl)
  DEALLOCATE (cvegh)
  DEALLOCATE (sncov)
  

!-------------------------------------------------------------------
! 3.  Read atmospheric temperature input  (Only CASE CITVEG = 'Tair' or CASE CIATM = 'Pellarin')
!-------------------------------------------------------------------


SELECT CASE (CITVEG)

   CASE ('Tair')

    OPEN(NULTMP,FILE='forcing_cmem_tair.asc',status='old')


     DO JJ = 1, N
      READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),ftair(JJ),doy(JJ)

     ENDDO
     CLOSE(NULTMP)

     ftveg(:) = ftair(:)

!    Update the mask to eliminate unrealistic points
     SELECT CASE (LOMASK_AUTO)
       CASE (.True.)
       WHERE (ftveg < 100.0_JPRM ) mask = 0_JPIM  ! Unrealistic Temperature
     END SELECT
!

   CASE ('Tsurf')

     ftveg(:) = ftl_lsm(:,1)

END SELECT

SELECT CASE (CIATM)

   CASE ('Pellarin', 'Ulaby')

    OPEN(NULTMP,FILE='forcing_cmem_tair.asc',status='old')


     DO JJ = 1, N
      READ (NULTMP,*) nbpt(JJ),ndate(JJ),ntime(JJ),ftair(JJ),doy(JJ)

     ENDDO
     CLOSE(NULTMP)

!    Update the mask to eliminate unrealistic points
      SELECT CASE (LOMASK_AUTO)
       CASE (.True.)
       WHERE (ftair < 100.0_JPRM ) mask = 0_JPIM  ! Unrealistic Temperature
      END SELECT
!


END SELECT

!-------------------------------------------------------------------
! 4.  Read SSS
!-------------------------------------------------------------------

sal_sea(:) = 32.5_JPRM

         
END SUBROUTINE RDCMEMASCII
