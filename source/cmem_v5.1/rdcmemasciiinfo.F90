! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE RDCMEMASCIIINFO

! Purpose :
! -------
!     Check consistency of ascii input data for cmem 
   
! Externals :
! ---------

! Internal datafields :
! -------------------
  
! Authors:
!     Patricia de Rosnay  01-10-2007
!     Patricia de Rosnay  19-01-2012 add HTessel confi with variable LAI
!------------------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRM
USE YOMLUN, ONLY : NULOUT, NULTMP
USE YOMCMEMPAR, ONLY : CITVEG, CIDVEG,CIATM, LGPRINT,LOFIELDEXP,nlay_soil_ls

USE YOMCMEMFIELDS, ONLY : N, JJ, NFCASCMAX

IMPLICIT NONE

CHARACTER(LEN=100) :: CHEAD
INTEGER(KIND=JPIM) :: i,inbpt,IC,idate,itime
REAL(KIND=JPRM) ::  zdoy,zfwc_lsm(nlay_soil_ls),ztskin,zftl_lsm(nlay_soil_ls),zfsnowd,zfrsnow &
                 &  ,zlaiL,zlaiH,ztveg,zsand,zclay &
                 &  , zfZ,zfTVL,zfTVH,zfvegl,zfvegh,zfwater


!----------------------------------------------------------
! 1.  Open main input file and check its size
!----------------------------------------------------------

OPEN(NULTMP,FILE='forcing_cmem_main.asc',status='old')
WRITE(NULOUT,*) 'Read forcing_cmem_main'


IC = 0_JPIM
DO JJ = 1,NFCASCMAX 

READ (NULTMP,*, END=300) inbpt,idate,itime,(zfwc_lsm(i),i=1,nlay_soil_ls),&
                         & ztskin,(zftl_lsm(i),i=1,nlay_soil_ls),zfsnowd,zfrsnow,zdoy
IC = IC + 1_JPIM

ENDDO 
300  CONTINUE
CLOSE(NULTMP)

N = IC
WRITE(NULOUT,*) 'Forcing file size, N= ', N

!----------------------------------------------------------
! 2.  Open soil and vegetation input files 
!----------------------------------------------------------


SELECT CASE (LOFIELDEXP)
       
     CASE (.TRUE.)
       ! constant parameters to be read
       ! no size check required

       OPEN(NULTMP,FILE='forcing_cmem_soil_constant.asc',status='old')
       WRITE(NULOUT,*) 'Read forcing_cmem_soil'
       READ (NULTMP,*) zsand,zclay, zfZ
       CLOSE(NULTMP)

       SELECT CASE (CIDVEG)

          CASE ('Tessel','HTessel')

          OPEN(NULTMP,FILE='forcing_cmem_veg-tessel_constant.asc',status='old')
          WRITE(NULOUT,*) 'Read forcing_cmem_veg'
          READ (NULTMP,*) zfTVL,zfTVH,zfvegl,zfvegh,zfwater
          CLOSE(NULTMP)
          
          CASE ('Ecoclimap')

          OPEN(NULTMP,FILE='forcing_cmem_veg-ecoclimap_constant.asc',status='old')
          WRITE(NULOUT,*) 'Read forcing_cmem_veg'
          READ (NULTMP,*) zfTVL,zfTVH,zfvegl,zfvegh,zfwater
          CLOSE(NULTMP)

        END SELECT

     CASE (.FALSE.)

       ! Surface conditions vary
        OPEN(NULTMP,FILE='forcing_cmem_soil.asc',status='old')
        WRITE(NULOUT,*) 'Read forcing_cmem_soil'

        IC = 0_JPIM
        DO JJ = 1,NFCASCMAX

         READ (NULTMP,*, END=301) inbpt,idate,itime,zsand,zclay,zfZ,zdoy
         IC = IC + 1_JPIM

        ENDDO
301     CONTINUE
        CLOSE(NULTMP)


        IF( IC /= N) CALL ABOR1('forcing_cmem_main.asc and forcing_cmem_soil.asc not consistent')
  


       SELECT CASE (CIDVEG)

         CASE ('Tessel','HTessel')

           OPEN(NULTMP,FILE='forcing_cmem_veg-tessel.asc',status='old')
           WRITE(NULOUT,*) 'Read forcing_cmem_veg-tessel'
           IC = 0_JPIM
           DO JJ = 1,NFCASCMAX
             READ (NULTMP,*, END=302) inbpt,idate,itime,zfTVL,zfTVH,zfvegl,zfvegh,zfwater,zdoy
             IC = IC + 1_JPIM
           ENDDO
302        CONTINUE
           CLOSE(NULTMP)
           IF( IC /= N) CALL ABOR1('forcing_cmem_main.asc and forcing_cmem_veg-tessel.asc not consistent')

         CASE ('Ecoclimap')

           OPEN(NULTMP,FILE='forcing_cmem_veg-ecoclimap.asc',status='old')
           WRITE(NULOUT,*) 'Read forcing_cmem_veg-ecoclimap'
           IC = 0_JPIM
           DO JJ = 1,NFCASCMAX
             READ (NULTMP,*, END=303) inbpt,idate,itime,zfTVL,zfTVH,zfvegl,zfvegh,zfwater,zdoy
             IC = IC + 1_JPIM
           ENDDO
303        CONTINUE
           CLOSE(NULTMP)
           IF( IC /= N) CALL ABOR1('forcing_cmem_main.asc and forcing_cmem_veg-ecoclimap.asc not consistent')

           CLOSE(NULTMP)

       END SELECT



  
END SELECT

!----------------------------------------------------------
! 3.  Open vegetation lai input file and check its size
!----------------------------------------------------------


SELECT CASE (CIDVEG)

   CASE ('Ecoclimap','HTessel')

    OPEN(NULTMP,FILE='forcing_cmem_lai.asc',status='old')
    WRITE(NULOUT,*) 'Read forcing_cmem_lai'


    IC = 0_JPIM
    DO JJ = 1,NFCASCMAX 

      READ (NULTMP,*, END=304) inbpt,idate,itime,zlaiL,zlaiH,zdoy
      IC = IC + 1_JPIM

    ENDDO
304 CONTINUE
    CLOSE(NULTMP)

IF( IC /= N) CALL ABOR1('Forcing files forcing_cmem_main.asc and forcing_cmem_lai-ecoclimap.asc not consistent')

   CASE ('Tessel')

   ! no input file required

END SELECT 



!--------------------------------------------------------------------------------------------------
! 4.  Open atmospheric temperature input file and check its size  (Only CASE CITVEG = 'Tair')
!--------------------------------------------------------------------------------------------------


SELECT CASE (CITVEG)

   CASE ('Tair')

   OPEN(NULTMP,FILE='forcing_cmem_tair.asc',status='old')
   WRITE(NULOUT,*) 'Read forcing_cmem_tair'

    IC = 0_JPIM
    DO JJ = 1,NFCASCMAX

     READ (NULTMP,*, END=307) inbpt,idate,itime,ztveg,zdoy
     IC = IC + 1_JPIM

    ENDDO
307 CONTINUE
    CLOSE(NULTMP)

    IF( IC /= N) CALL ABOR1('Forcing forcing_cmem_main.asc and forcing_cmem_tair.asc not consistent')
END SELECT
!
SELECT CASE (CIATM)

   CASE ('Pellarin', 'Ulaby')
   OPEN(NULTMP,FILE='forcing_cmem_tair.asc',status='old')
   WRITE(NULOUT,*) 'Read forcing_cmem_tair'

    IC = 0_JPIM
    DO JJ = 1,NFCASCMAX

     READ (NULTMP,*, END=308) inbpt,idate,itime,ztveg,zdoy
     IC = IC + 1_JPIM

    ENDDO
308 CONTINUE
    CLOSE(NULTMP)

    IF( IC /= N) CALL ABOR1('Forcing forcing_cmem_main.asc and forcing_cmem_tair.asc not consistent')

END SELECT





WRITE(NULOUT,*) 'Forcing_info: input files are consistent. Grid dimension: ',N
         
END SUBROUTINE RDCMEMASCIIINFO
