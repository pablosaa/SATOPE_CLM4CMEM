! Copyright 2006-2014 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities granted to it by
! virtue of its status as an intergovernmental organisation nor does it submit to any jurisdiction.
!
SUBROUTINE WRCMEMGRIBAPI

! Purpose :
! -------
!  WRITE CMEM OUTPUT DATA IN NETCDF FORMAT
   
! Interface :
! ---------

! Method :
! ------

! Externals :
! ---------

! Internal variables
! ------------------

! Author :
!     Patricia de Rosnay  January 2009
!     PdR: May 2012: add angle info in output file name

!------------------------------------------------------------------------------
!
USE grib_api
!
USE PARKIND1, ONLY : JPIM, JPRM 
USE YOMLUN, ONLY : NULOUT,NULTMP
!
USE YOMCMEMPAR, ONLY : JPHISTLEV,CNAMEID,CFREQ,CANGLE
USE YOMCMEMFIELDS, ONLY: N,JJ,CLNAME,fwc_veg,fb,ftfrac,fh,fWP,ftau_atm,ftb_au &
                     & ,ftb_toa,fteffC,fsurf_emis,ftau_veg
USE YOMCMEMGRIBAPI
!
!
IMPLICIT NONE
!
! Local variables:
!
INTEGER(KIND=JPIM) :: igribtemplate,igribclone,igribinput,outfile
INTEGER(KIND=JPIM) ::  ifile,iret,date,indicatorOfParameter

REAL(KIND=JPRM), ALLOCATABLE :: field2D(:)
REAL(KIND=JPRM)              :: missingValue

missingValue = 9999.0_JPRM


!
! 1.0 Template and reference
!---------------------------
!
! 1.1 Get reference file
!-----------------------
!
! Open the verification file

CALL grib_open_file(ifile,'SWVL1.grib','r')

! Load messages

CALL grib_new_from_file(ifile,igribinput, iret)

! Close the verification input file

CALL grib_close_file(ifile)

CALL grib_get_size(igribinput,'values',numberOfPoints)

!
! 1.2 Get template file
!----------------------
!
CALL grib_new_from_template(igribtemplate, "SWVL1")

CALL grib_get_size(igribtemplate,'values',numberOfValues)

WRITE(NULOUT,*) 'numberOfValues=',numberOfValues
IF( ALLOCATED (field2D)) DEALLOCATE (field2D)
ALLOCATE(field2D(numberOfValues), stat=iret)

CALL grib_get(igribtemplate,'dataDate',date)
WRITE(NULOUT,*) "Date = ",date


!
! 1.3 Check consistency between template file and input file
!-----------------------------------------------------------
!

IF (numberOfValues /= numberOfPoints) CALL ABOR1('WRCMEMGRIBAPI: Template size not match')
    
! Release data (frees up any internal GRIBAPI resources)

CALL grib_release(igribinput)


! 2. WRITE outputs
!---------------------------------------
!
!

HISTLEV: SELECT CASE (JPHISTLEV)

! Level 1 outputs:
!-----------------

  CASE ( 1_JPIM, 2_JPIM, 3_JPIM )

!
!  TBH
!-------
 !CLNAME='out_tbtoaH_'//CNAMEID//'.grib'
  CLNAME='out_tbtoaH_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 194 
    field2D(:) =  ftb_toa(:,1)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)
!
!
!  TBV
!--------
! CLNAME='out_tbtoaV_'//CNAMEID//'.grib'
  CLNAME='out_tbtoaV_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 194 
    field2D(:) = ftb_toa(:,2)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
!
!  TEFF
!-------
! CLNAME='out_tmpeff_'//CNAMEID//'.grib'
  CLNAME='out_tmpeff_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 130 
    field2D(:) = fteffC(:,2)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)




  SELECT CASE (JPHISTLEV)

! Level 2 outputs:
!-----------------

  CASE ( 2_JPIM, 3_JPIM )

!
!  tau veg H
!-----------
! CLNAME='out_tauvgH_'//CNAMEID//'.grib'
  CLNAME='out_tauvgH_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) = ftau_veg(:,1)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
!  tau veg V
!-----------
! CLNAME='out_tauvgV_'//CNAMEID//'.grib'
  CLNAME='out_tauvgV_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) = ftau_veg(:,2)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
! Bare fract 
!-----------
! CLNAME='out_frbare_'//CNAMEID//'.grib'
  CLNAME='out_frbare_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) = ftfrac(:,1)+ftfrac(:,2)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
! VWC sum
!-----------
! CLNAME='out_VWCsum_'//CNAMEID//'.grib'
  CLNAME='out_VWCsum_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) = fwc_veg(:,1)*(ftfrac(:,3)+ftfrac(:,4)) +  fwc_veg(:,2)*(ftfrac(:,5)+ftfrac(:,6))


! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
! tau atmosphere
!---------------
! CLNAME='out_tau_at_'//CNAMEID//'.grib'
  CLNAME='out_tau_at_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) = ftau_atm(:)


! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
! TB atm upward 
!---------------
! CLNAME='out_tb_aup_'//CNAMEID//'.grib'
  CLNAME='out_tb_aup_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) = ftb_au(:)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)



  SELECT CASE (JPHISTLEV)


! LEVEL 3 outputs:
!------------------
  CASE ( 3_JPIM )

!
! Bare soil fraction (without snow)
!----------------------------------
! CLNAME='out_ffrac1_'//CNAMEID//'.grib'
  CLNAME='out_ffrac1_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) = ftfrac(:,1)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
! Low veg fraction (without snow)
!--------------------------------
! CLNAME='out_ffrac3_'//CNAMEID//'.grib'
  CLNAME='out_ffrac3_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) = ftfrac(:,3)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
! VWC of Low Vegetation 
!----------------------
! CLNAME='out_VWC_LO_'//CNAMEID//'.grib'
  CLNAME='out_VWC_LO_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) = fwc_veg(:,1)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
! VWC of High Vegetation 
!------------------------
! CLNAME='out_VWC_HI_'//CNAMEID//'.grib'
  CLNAME='out_VWC_HI_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) = fwc_veg(:,2)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
! B parameter for low veg
!------------------------
! CLNAME='out_Bpa_LO_'//CNAMEID//'.grib'
  CLNAME='out_Bpa_LO_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) = fb(:,1)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
! B parameter for high veg
!-------------------------
! CLNAME='out_Bpa_HI_'//CNAMEID//'.grib'
  CLNAME='out_Bpa_HI_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) = fb(:,2)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
! Roughness param h
!------------------
! CLNAME='out_fh_par_'//CNAMEID//'.grib'
  CLNAME='out_fh_par_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) = fh(:)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
! Wilting point soil moisture
!----------------------------
! CLNAME='out_WP_par_'//CNAMEID//'.grib'
  CLNAME='out_WP_par_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
   field2D(:) =  fWP(:)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
! Emissivity at H pol 
!--------------------
! CLNAME='out_emissH_'//CNAMEID//'.grib'
  CLNAME='out_emissH_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 124 
    field2D(:) =  fsurf_emis(:,1)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

!
! Emissivity at V pol 
!--------------------
! CLNAME='out_emissV_'//CNAMEID//'.grib'
  CLNAME='out_emissV_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 124 
    field2D(:) =  fsurf_emis(:,2)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)

! C paramareter for Teff
!-----------------------
! CLNAME='out_teff_C_'//CNAMEID//'.grib'
  CLNAME='out_teff_C_'//CNAMEID//'_'//cangle//'.grib'

! ALLOCATE and get values
    indicatorOfParameter = 80 
    field2D(:) =  fteffC(:,1)

! Open clone the template and define output file based on cloned template
    CALL grib_open_file(outfile, CLNAME,'w')
    CALL grib_clone(igribtemplate,igribclone)
    CALL grib_set(igribclone,"dataDate", date)
    CALL grib_set(igribclone,"indicatorOfParameter", indicatorOfParameter)
    CALL grib_set(igribclone, 'missingValue',missingValue)

! Write output file and close
    CALL grib_set(igribclone,'values',field2D)
    CALL grib_write(igribclone,outfile)
    CALL grib_close_file(outfile)
! Release template
    CALL grib_release(igribclone)


ENDSELECT
ENDSELECT
ENDSELECT HISTLEV



! Release template
CALL grib_release(igribtemplate)
! Clean field2D
DEALLOCATE(field2D)


       
END SUBROUTINE WRCMEMGRIBAPI

