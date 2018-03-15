! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
MODULE YOMCMEMGRIBAPI

! Module containing GRIB API file handlers

!---------------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,JPRM

IMPLICIT NONE


REAL(KIND = JPRM), ALLOCATABLE   ::  fapifield(:,:)
REAL(KIND = JPRM)                ::  latitudeOfFirstgridPointInDegrees
REAL(KIND = JPRM)                ::  longitudeOfFirstgridPointInDegrees
REAL(KIND = JPRM)                ::  latitudeOfLastgridPointInDegrees
REAL(KIND = JPRM)                ::  longitudeOfLastgridPointInDegrees
INTEGER(KIND=JPIM)               ::  numberOfPointsAlongAParallel
INTEGER(KIND=JPIM)               ::  numberOfPointsAlongAMeridian
INTEGER(KIND=JPIM)               ::  numberOfLevels
INTEGER(KIND=JPIM)               ::  numberOfPoints,numberOfValues
CHARACTER(LEN=10)  CCVARAPI_NAME






END MODULE YOMCMEMGRIBAPI
