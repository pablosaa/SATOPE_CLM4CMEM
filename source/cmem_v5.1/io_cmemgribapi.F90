! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE IO_CMEMGRIBAPI(CLNAME,LD_firstcall)

! PURPOSE
! -------
!  To read data from / to GRIBAPI files
   
! INTERFACE
! ---------
!  CALL IO_CMEMGRIBAPI

! METHOD
! ------
!
! EXTERNALS
! ---------
!
! INTERNAL ROUTINES
! -----------------
!  Patricia de Rosnay  December 2008
!------------------------------------------------------------------------------

USE YOMCMEMGRIBAPI

IMPLICIT NONE

CHARACTER(len=80) :: clname
LOGICAL :: LD_firstcall
!------------------------------------------------------------------------------


    CALL INGRIBAPI(CLNAME,LD_firstcall)

RETURN

CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE INGRIBAPI(INFILE,LD_firstcall)

USE grib_api
USE YOMLUN, ONLY : NULOUT
USE YOMCMEMGRIBAPI
USE YOMCMEMPAR, ONLY : LGPRINT
!
IMPLICIT NONE


CHARACTER(len=80) infile

LOGICAL :: LD_firstcall
INTEGER(KIND=JPIM) ::  ifile,iret,i
REAL, ALLOCATABLE :: values(:)

! local grib file (dimension on number of levels):
INTEGER,DIMENSION(:),ALLOCATABLE   ::  igrib



! 1.  Open the input file
!------------------------

 WRITE(NULOUT,*)'INPUT_GRIBAPI, read file  ',trim(infile)
 call grib_open_file(ifile,INFILE,'r') 



!
! 2. Count number of levels, load message(s) from the grib file
!    and close the file
!----------------------------------------------------------------
!

! count levels
  call grib_count_in_file(ifile,numberOfLevels)
  IF( ALLOCATED (igrib)) DEALLOCATE (igrib)
  ALLOCATE(igrib(numberOfLevels))
  igrib=-1

! load messages
  DO i=1,numberOfLevels
     call grib_new_from_file(ifile,igrib(i), iret)
  END DO

 ! close the file
  call grib_close_file(ifile)


!
! 3. Get dimension of the field GRIBAPI variables
!----------------------------------------------------------------
!
 ! Loop on all the messages in memory
  DO i=1,numberOfLevels


     call grib_get(igrib(i),'numberOfPointsAlongAMeridian', &
          numberOfPointsAlongAMeridian)
!     write(*,*) 'numberOfPointsAlongAMeridian=', &
!          numberOfPointsAlongAMeridian

     call grib_get(igrib(i), 'latitudeOfFirstGridPointInDegrees', &
          latitudeOfFirstgridPointInDegrees)
!     write(*,*) 'latitudeOfFirstGridPointInDegrees=', &
!          latitudeOfFirstgridPointInDegrees

     call grib_get(igrib(i), 'longitudeOfFirstGridPointInDegrees', &
          longitudeOfFirstgridPointInDegrees)
!     write(*,*) 'longitudeOfFirstGridPointInDegrees=', &
!          longitudeOfFirstgridPointInDegrees

     call grib_get(igrib(i), 'latitudeOfLastGridPointInDegrees', &
          latitudeOfLastgridPointInDegrees)
!     write(*,*) 'latitudeOfLastGridPointInDegrees=', &
!          latitudeOfLastgridPointInDegrees

     call grib_get(igrib(i), 'longitudeOfLastGridPointInDegrees', &
          longitudeOfLastgridPointInDegrees)

    !     get the size of the values array
     call grib_get_size(igrib(i),'values',numberOfPoints)
    ! write(*,*) 'numberOfPoints=',numberOfPoints


ENDDO 

!
!
! 3. Allocate and Read input variable according to its dimensions 
!-----------------------------------------------------------------
!
SELECT CASE (LD_firstcall)
!
  CASE(.False.)
!
   IF( ALLOCATED (fapifield)) DEALLOCATE (fapifield)
   ALLOCATE(fapifield(numberOfPoints,numberOfLevels))
!
!
  DO i=1,numberOfLevels
    IF( ALLOCATED (values)) DEALLOCATE (values)
    ALLOCATE(values(numberOfPoints), stat=iret)
     !     get data values
     call grib_get(igrib(i),'values',values)

     fapifield(:,i) = values(:) 
  ENDDO
!
END SELECT

!
!
! 5. Deallocate local variables and close the file. This frees up any internal GRIBAPI resources
!----------------------------------------------------------------------------------------------
!

 DO i=1,numberOfLevels
    call grib_release(igrib(i))
  END DO

IF( ALLOCATED (values)) DEALLOCATE(values)
IF( ALLOCATED (igrib)) DEALLOCATE(igrib)


RETURN


END SUBROUTINE INGRIBAPI

!------------------------------------------------------------------------------


END SUBROUTINE IO_CMEMGRIBAPI


