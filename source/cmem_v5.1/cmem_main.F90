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


!!$
!!$CALL CMEM_SETUP
!!$!
!!$!
!!$! 1. Get information on the input files and allocate variables accordingly
!!$!------------------------------------------------------------------------------
!!$!
!!$! Read input file and get information on grid size N
!!$!
!!$!
!!$SELECT CASE ( CFINOUT )
!!$!
!!$  CASE ( 'gribapi' ) !gribapicase
!!$!
!!$! PSG commented.   CALL RDCMEMGRIBAPIINFO !gribapicase
!!$!
!!$  CASE ( 'ascii' )
!!$!
!!$  CALL RDCMEMASCIIINFO
!!$!
!!$  CASE ( 'netcdf' ) !netcdfcase
!!$!
!!$  CALL RDCMEMNETCDFINFO !netcdfcase
!!$!
!!$  CASE ('clm')  ! PSG: new case for CLM support     
!!$     call info_CLM_file(CLNAME,Ntot=Nvars,inhr=1)
!!$     NLONS = Nvars(1)
!!$     NLATS = Nvars(2)
!!$     NTIMES = Nvars(3)
!!$     N = product(Nvars)
!!$!
!!$!  CASE ( 'ifs' )
!!$!!
!!$!  N = 1
!!$!
!!$END SELECT
!!$!
!!$!
!!$!
!!$! 1.2 Allocate variables according to the size of input
!!$!
!!$ALLOCATE (ftau_atm(N))
!!$ALLOCATE (ftau_veg(N,2))
!!$ALLOCATE (ftb_au(N))
!!$ALLOCATE (ftb_ad(N))
!!$ALLOCATE (ftair(N))
!!$ALLOCATE (fsnowd(N))
!!$ALLOCATE (frsnow(N))
!!$ALLOCATE (fwater(N))
!!$ALLOCATE (mask(N))
!!$ALLOCATE (fTVL(N))
!!$ALLOCATE (fTVH(N))
!!$ALLOCATE (fs_laiL(N))
!!$ALLOCATE (ftfrac(N,JPCMEMTILE))
!!$ALLOCATE (ftskin(N))
!!$ALLOCATE (ftl_lsm(N,nlay_soil_ls))
!!$ALLOCATE (fwc_lsm(N,nlay_soil_ls))
!!$ALLOCATE (fsand(N))
!!$ALLOCATE (fclay(N))
!!$ALLOCATE (frho_b(N))
!!$ALLOCATE (fp(N))
!!$ALLOCATE (fWP(N))
!!$ALLOCATE (falpha(N))
!!$ALLOCATE (fsal_wat(N))
!!$ALLOCATE (fh(N))
!!$ALLOCATE (ftveg(N))
!!$ALLOCATE (fwc_veg(N,2))
!!$ALLOCATE (fb(N,2))
!!$ALLOCATE (ftauN(N,2))
!!$ALLOCATE (fw_effL(N,2))
!!$ALLOCATE (fw_effH(N,2))
!!$ALLOCATE (tsoil(nlay_soil_mw))
!!$ALLOCATE (fNrh_L(N))
!!$ALLOCATE (fNrh_H(N))
!!$ALLOCATE (fNrv_L(N))
!!$ALLOCATE (fNrv_H(N))
!!$ALLOCATE(fhrmodel(N,2))
!!$ALLOCATE (ftth(N,2))
!!$ALLOCATE (fttv(N,2))
!!$ALLOCATE (fZ(N))
!!$ALLOCATE (sal_sea(N))
!!$ALLOCATE(ftheta_inc(N,NINC))  ! PSG: allocation of pixel-based theta
!!$ALLOCATE(TB_HV(NLONS,NLATS,2,NINC,NTIMES)) ! PSG: TB matrix for sat_OP
!!$print*, 'NINC is: ', NINC
!!$!
!!$!
!!$SELECT CASE (CFINOUT)
!!$   CASE ('ascii')
!!$   ALLOCATE (ndate(N))
!!$   ALLOCATE (ntime(N))
!!$   ALLOCATE (nbpt(N))
!!$   ALLOCATE (doy(N))
!!$END SELECT
!!$!
!!$!  allocate output
!!$!  ----------------
!!$!
!!$ALLOCATE ( ftb_soil(N,2) )
!!$ALLOCATE ( ftb_toa(N,2) )
!!$ALLOCATE ( fsurf_emis(N,2) )
!!$ALLOCATE ( fteffC(N,3) )
!!$!
!!$!
!!$! 2. Read input data
!!$!      ---------------
!!$!
!!$SELECT CASE (CFINOUT)
!!$!
!!$CASE ('gribapi') !gribapicase
!!$!
!!$   print*,'not supported' ! PSG commented.   CALL RDCMEMGRIBAPI !gribapicase
!!$!
!!$CASE ('ascii')
!!$   !
!!$   CALL RDCMEMASCII
!!$   !
!!$CASE ('netcdf') !netcdfcase
!!$!
!!$   CALL RDCMEMNETCDF !netcdfcase
!!$!
!!$CASE ('clm')   ! PSG: following 3 lines clm4cmem implementations:
!!$   call read_CLM_file(CLNAME,inhr=1,SAT=SAT,LS=CLMVARS)
!!$   call memory_cmem_forcing(CLMVARS)
!!$
!!$!  CASE ('ifs')
!!$!!
!!$!   CALL RDCMEMIFS
!!$!
!!$END SELECT
!!$!
!!$!
!!$WRITE(NULOUT,*) 'Read input files done'
!!$!
!!$IANGLE: DO JJINC = 1_JPIM,NINC  ! PSG: here start loop over incidence angles
!!$! 2.1 Compute atmospheric and vegetation optical depth,  and soil roughness
!!$!
!!$CALL CMEM_INIT
!!$!
!!$!
!!$! 2.0 deallocate variables that are not useful anymore
!!$!
!!$DEALLOCATE (fwater) ! read in RD_frac used in RD_soil
!!$DEALLOCATE (fTVL)   ! read in RD_frac used in RD_veg
!!$DEALLOCATE (fTVH)   ! read in RD_frac used in RD_veg
!!$DEALLOCATE (fs_laiL) ! read in RD_frac used in RD_veg
!
!!$DEALLOCATE (fZ)
!!$!
!!$!
!!$!
!!$IF (LGPRINT) WRITE(NULOUT,*) 'CMEM_main, end of 1.2'
!!$!
!!$!
!!$! 3. Compute toa brightness temperatures
!!$!    -----------------------------------
!!$!
!!$WRITE(NULOUT,*) 'CMEM_main, compute TB field ', SAT%theta(JJINC) ! PSG:
!!$!
!!$!
!!$
!!$FIELD: DO JJ = 1, N
!!$!
!!$! 3.1 Initialize
!!$!---------------
!!$!
!!$tb_soil(:) = (/0.,0./)
!!$tb_veg(:) = (/0.,0./)
!!$tb_toa(:) = (/0.,0./)
!!$fsurf_emis(JJ,:) = (/0.,0./)
!!$ftb_soil(JJ,:) = (/0.,0./)
!!$ftau_veg(JJ,:) = (/0.,0./)
!!$costheta = cos(pi*ftheta_inc(JJ,JJINC)/180._JPRM)  ! PSG: pixel-wise incidence
!!$sintheta = sin(pi*ftheta_inc(JJ,JJINC)/180._JPRM)  ! PSG: pixel-wise incidence
!!$!
!!$SELECT CASE (MASK(JJ))
!!$!
!!$  CASE ( 0_JPIM )    ! points on which TB is not computed (merge MASK_AUTO and MASK_OCEAN)
!!$!
!!$    CYCLE
!!$!
!!$  CASE ( 1_JPIM )    ! points on which TB is computed
!!$!
!!$! Cell values
!!$  tfrac(:) = (ftfrac(JJ,:))
!!$  tau_atm = ftau_atm(JJ) / costheta
!!$  tb_au = ftb_au(JJ)
!!$  tb_ad = ftb_ad(JJ)
!!$  IF (LGPRINT)  WRITE(NULOUT,*) '--- TBSKY up:' ,tb_au
!!$  IF (LGPRINT)  WRITE(NULOUT,*) '--- TBSKY down:' ,tb_ad
!!$  t_veg = ftveg(JJ)
!!$  tair = ftair(JJ)
!!$!
!!$!
!!$  DO JCMEMTILE = 1, JPCMEMTILE
!!$!
!!$  IF (LGPRINT) WRITE(NULOUT,*) '--- Mean TEFF a',fteffC(JJ,:)
!!$    IF (ftfrac(JJ,JCMEMTILE) == 0.)   CYCLE
!!$!
!!$IF (LGPRINT) WRITE(NULOUT,*) '----------- Tile',JCMEMTILE
!!$!
!!$!
!!$! 3.2 Compute surface emissivity
!!$!
!!$    SELECT CASE (JCMEMTILE)
!!$      CASE ( 1,2 )
!!$        Nrh = fNrh_L(JJ)
!!$        Nrv = fNrv_L(JJ)
!!$        hrmodel = fhrmodel(JJ,1)
!!$        tsoildeep = ftl_lsm(JJ,nlay_soil_ls)
!!$        CALL CMEM_SOIL
!!$        IF (LGPRINT) WRITE(NULOUT,*) '--- Mean TEFF bf',&
!!$          & ftfrac(JJ,JCMEMTILE),fteffC(JJ,1),t_eff(1),fteffC(JJ,2),t_eff(2),&
!!$          & fteffC(JJ,3),t_eff(3)
!!$        fteffC(JJ,:) = fteffC(JJ,:) + ftfrac(JJ,JCMEMTILE) * t_eff(:)
!!$        IF (LGPRINT) WRITE(NULOUT,*) '--- Mean TEFF b',fteffC(JJ,:)
!!$      CASE ( 6 )
!!$        Nrh = fNrh_H(JJ)
!!$        Nrv = fNrv_H(JJ)
!!$        hrmodel = fhrmodel(JJ,2)
!!$        tsoildeep = ftl_lsm(JJ,1)
!!$        CALL CMEM_SOIL
!!$        IF (LGPRINT) WRITE(NULOUT,*) '--- Mean TEFF bf',&
!!$          & ftfrac(JJ,JCMEMTILE),fteffC(JJ,1),t_eff(1),fteffC(JJ,2),t_eff(2),&
!!$          & fteffC(JJ,3),t_eff(3)
!!$        fteffC(JJ,:) = fteffC(JJ,:) + ftfrac(JJ,JCMEMTILE) * t_eff(:)
!!$        IF (LGPRINT) WRITE(NULOUT,*) '--- Mean TEFF c',fteffC(JJ,:)
!!$      CASE (  3, 4 )
!!$        Nrh = fNrh_L(JJ)
!!$        Nrv = fNrv_L(JJ)
!!$        hrmodel = fhrmodel(JJ,1)
!!$        tsoildeep = ftl_lsm(JJ,nlay_soil_ls)
!!$        CALL CMEM_SOIL
!!$        IF (LGPRINT) WRITE(NULOUT,*) '--- Mean TEFF bf',&
!!$          & ftfrac(JJ,JCMEMTILE),fteffC(JJ,1),t_eff(1),fteffC(JJ,2),t_eff(2),&
!!$          & fteffC(JJ,3),t_eff(3)
!!$        fteffC(JJ,:) = fteffC(JJ,:) + ftfrac(JJ,JCMEMTILE) * t_eff(:)
!!$        IF (LGPRINT) WRITE(NULOUT,*) '--- Mean TEFF d',fteffC(JJ,:)
!!$      CASE ( 5 )
!!$        Nrh = fNrh_H(JJ)
!!$        Nrv = fNrv_H(JJ)
!!$        hrmodel = fhrmodel(JJ,2)
!!$        tsoildeep = ftl_lsm(JJ,nlay_soil_ls)
!!$        CALL CMEM_SOIL
!!$        IF (LGPRINT) WRITE(NULOUT,*) '--- Mean TEFF bf',&
!!$          & ftfrac(JJ,JCMEMTILE),fteffC(JJ,1),t_eff(1),fteffC(JJ,2),t_eff(2),&
!!$          & fteffC(JJ,3),t_eff(3)
!!$        fteffC(JJ,:) = fteffC(JJ,:) + ftfrac(JJ,JCMEMTILE) * t_eff(:)
!!$        IF (LGPRINT) WRITE(NULOUT,*) '--- Mean TEFF e',fteffC(JJ,:)
!!$      CASE ( 7 )
!!$        Nrh = 0.0_JPRM
!!$        Nrv = 0.0_JPRM
!!$        hrmodel = 0.0_JPRM
!!$        tsoildeep = ftl_lsm(JJ,nlay_soil_ls)
!!$        CALL CMEM_SOIL
!!$        IF (LGPRINT) WRITE(NULOUT,*) '--- Mean TEFF bf',&
!!$          & ftfrac(JJ,JCMEMTILE),fteffC(JJ,1),t_eff(1),fteffC(JJ,2),t_eff(2),&
!!$          & fteffC(JJ,3),t_eff(3)
!!$        !fteffC(JJ,1) = fteffC(JJ,1) + ftfrac(JJ,JCMEMTILE) * t_eff(1)
!!$        fteffC(JJ,2) = fteffC(JJ,2) + ftfrac(JJ,JCMEMTILE) * t_eff(2)
!!$        fteffC(JJ,3) = fteffC(JJ,3) + ftfrac(JJ,JCMEMTILE) * t_eff(3)
!!$        IF (LGPRINT) WRITE(NULOUT,*) '--- Mean TEFF f',fteffC(JJ,:)
!!$    END SELECT
!!$!
!!$! 3.3 Compute vegetation emissivity if present
!!$!
!!$!
!!$    SELECT CASE (JCMEMTILE)
!!$      CASE ( 3, 4 )  ! Low vegetation
!!$        w_eff(:) = fw_effL(JJ,:)
!!$        a_geo = a_geoL
!!$        wc_veg =  fwc_veg(JJ,1)
!!$        bj =  fb(JJ,1)
!!$        tauN = ftauN(JJ,1)
!!$        tth = ftth(JJ,1)
!!$        ttv = fttv(JJ,1)
!!$        IF (tb_veg(1) == 0.)  CALL CMEM_VEG
!!$     CASE ( 5 )  ! High vegetation
!!$        w_eff(:) = fw_effH(JJ,:)
!!$        a_geo = a_geoH
!!$        wc_veg = fwc_veg(JJ,2)
!!$        bj =  fb(JJ,2)
!!$        tauN = ftauN(JJ,2)
!!$        tth = ftth(JJ,2)
!!$        ttv = fttv(JJ,2)
!!$        CALL CMEM_VEG
!!$     CASE ( 1, 2, 6, 7 )
!!$        ! high veg (case 6) will be accounted for on the top of snow in 3.6
!!$        ! so, for now do not account for high veg vegetation
!!$        tau_veg(:) = (/0.,0./)
!!$        tb_veg(:) = (/0.,0./)
!!$    END SELECT
!!$!
!!$!
!!$! 3.4 Compute top-of-vegetation brightness temperature
!!$!
!!$    CALL CMEM_RTM
!!$!
!!$!  3.5 Add snow layer: tb_tov-->tb top-of-snow
!!$!
!!$    SELECT CASE (JCMEMTILE)
!!$      CASE ( 2_JPIM, 4_JPIM, 6_JPIM )
!!$        XMV = 0.1 ! Snow Moisture, should be from data field
!!$        DO JJPOL = 1,2    ! h- and v-polarization
!!$          RSN(JJPOL)=1._JPRM  - tb_tov(JJPOL) / ftl_lsm(JJ,1)
!!$        ENDDO
!!$        CALL CMEM_SNOW(RSN,frsnow(JJ)/1000.,fsnowd(JJ),ESN)
!!$    END SELECT
!!$!
!!$!
!!$! 3.6 Add vegetation layer for High vegetation over snow
!!$!
!!$    SELECT CASE (JCMEMTILE)
!!$      CASE ( 6_JPIM )  ! High vegetation, snow
!!$        ! soil
!!$        DO JJPOL = 1,2    ! h- and v-polarization
!!$          r_r(JJPOL)=1. - ESN(JJPOL)
!!$        ENDDO
!!$        tb_soil = tb_tov
!!$        ! vegetation
!!$        w_eff(:) = fw_effH(JJ,:)
!!$        a_geo = a_geoH
!!$        wc_veg = fwc_veg(JJ,2)
!!$        bj =  fb(JJ,2)
!!$        CALL CMEM_VEG
!!$        ! Tb top-of-vegetation
!!$        CALL CMEM_RTM
!!$    END SELECT
!!$!
!!$! 3.7 Compute contribution to top-of-atmosphere brightness temperature
!!$!
!!$     tb_toa(:) = tb_toa(:) + ftfrac(JJ,JCMEMTILE) * tb_tov(:)
!!$     fsurf_emis(JJ,:) = fsurf_emis(JJ,:) + ftfrac(JJ,JCMEMTILE) * surf_emis(:)
!!$     ftb_soil(JJ,:) = ftb_soil(JJ,:) + ftfrac(JJ,JCMEMTILE) * tb_soil(:)
!!$!    Diagnostic Tau_veg
!!$     ftau_veg(JJ,:) = ftau_veg(JJ,:) + ftfrac(JJ,JCMEMTILE) * tau_veg(:)
!!$!    Diagnostic mean TEFF
!!$     IF (LGPRINT) WRITE(NULOUT,*) 'end of cmem tile loop',JCMEMTILE
!!$     IF (LGPRINT) WRITE(NULOUT,*) '                 TB TOV',tb_tov(:)
!!$     IF (LGPRINT) WRITE(NULOUT,*) '                 tb_soil',tb_soil(:)
!!$
!!$!
!!$  ENDDO ! end JCMEMTILE
!!$!
!!$!
!!$!   3.7 Top-of-Atmosphere brightness temperature
!!$!
!!$  IF (LGPRINT) WRITE(NULOUT,*) '--- TB no atm:',tb_toa(:)
!!$  ftb_toa(JJ,:) = tb_toa(:) * exp(-tau_atm) + ftb_au(JJ)
!!$  IF (LGPRINT) WRITE(NULOUT,*) '--- TBTOA:',ftb_toa(JJ,:)
!!$!
!!$!
!!$END SELECT  ! end of mask selection
!!$!
!!$ENDDO FIELD
!!$!
!!$! 3.8 PSG: In case CLM process the high-res results to produce SAT-res
!!$! --------------------------------------------------------------
!!$if(CFINOUT.eq.'clm'.and.allocated(TB_HV)) then
!!$   TB_HV(:,:,1,JJINC,:)=reshape(ftb_toa(:,1),(/NLONS,NLATS,NTIMES/))
!!$   TB_HV(:,:,2,JJINC,:)=reshape(ftb_toa(:,2),(/NLONS,NLATS,NTIMES/))
!!$   ! PSG: when last incidence angle, calculate satellite-like data
!!$   if(JJINC.eq.NINC) then
!!$      call TBSAT_OPERATOR(SAT,CLMVARS%lons,CLMVARS%lats,&
!!$           CLMVARS%resol_km,TB_HV)
!!$   end if
!!$end if
!!$WRITE(CANGLE,'(I2)') INT(SAT%theta(JJINC))  ! PSG: changing theta
!!$
!!$!
!!$! 4. Write outputs
!!$!-----------------
!!$!
!!$WRITE(NULOUT,*) 'CMEM_main, write output'
!!$!
!!$
!!$SELECT CASE (CFINOUT)
!!$!
!!$  CASE ('gribapi') !gribapicase
!!$!
!!$! PSG commented.     CALL WRCMEMGRIBAPI !gribapicase
!!$!
!!$  CASE ('ascii')
!!$!
!!$    CALL WRCMEMASCII
!!$!
!!$  CASE ('netcdf') !netcdfcase
!!$!
!!$     CALL WRCMEMNETCDF !netcdfcase
!!$!
!!$  CASE ('clm') ! PSG: CLM case
!!$     if(JPHISTLEV.lt.4_JPIM) then
!!$        ! PSG: Writing High-res NetCDF only.
!!$        CALL WRCMEMNETCDF
!!$     else if(JPHISTLEV.eq.4_JPIM) then
!!$        ! PSG: Writing Satellite Operator NetCDF only.
!!$        CLNAME='../output/out_level4_'//CNAMEID//'_'//cfreq//'_'//SAT%name//'.nc'
!!$        call write_satellite_operator(SAT,CLNAME)
!!$     else if(JPHISTLEV.eq.5_JPIM) then
!!$        ! PSG: Writing High-res level1 AND satellite opertor NetCDFs.
!!$        CALL WRCMEMNETCDF ! PSG: for netcdf
!!$        CLNAME='../output/out_level4_'//CNAMEID//'_'//cfreq//'_'//SAT%name//'.nc'
!!$        call write_satellite_operator(SAT,CLNAME)
!!$     else if(JPHISTLEV.eq.6_JPIM) then
!!$        WRITE(NULOUT,*) 'Data to keep in memory! no file storated!'
!!$     else
!!$        WRITE(NULOUT,*) 'ERROR by selecting the write JPHISTLEV param.'
!!$     end if
!!$
!!$
!!$!
!!$!  CASE ('ifs')
!!$!!
!!$!    CALL WRCMEMIFS
!!$!
!!$END SELECT
!!$!
!!$ENDDO IANGLE ! PSG: end over loop theta_inc
!!$
!!$!
!!$! 5 Clean up
!!$!   --------
!!$!
!!$! -- PSG: deallocation block coming from after CMEM_INIT
!!$DEALLOCATE (fwater) ! read in RD_frac used in RD_soil  PSG:
!!$DEALLOCATE (fTVL)   ! read in RD_frac used in RD_veg   PSG:
!!$DEALLOCATE (fTVH)   ! read in RD_frac used in RD_veg   PSG:
!!$DEALLOCATE (fs_laiL) ! read in RD_frac used in RD_veg  PSG:
!!$!
!!$DEALLOCATE (fZ)   ! PSG:
!!$! -- PSG: end of deallocations from location after CMEM_INIT
!!$DEALLOCATE ( fteffC )
!!$!
!!$DEALLOCATE (ftb_soil)
!!$DEALLOCATE (fsurf_emis)
!!$!
!!$DEALLOCATE (ftveg)
!!$DEALLOCATE (ftskin)
!!$DEALLOCATE (ftl_lsm)
!!$!
!!$DEALLOCATE (fwc_lsm)
!!$DEALLOCATE (ftb_toa)
!!$!
!!$DEALLOCATE (fsand)
!!$DEALLOCATE (fclay)
!!$DEALLOCATE (mask)
!!$DEALLOCATE (frho_b)
!!$DEALLOCATE (fp)
!!$DEALLOCATE (fWP)
!!$DEALLOCATE (falpha)
!!$DEALLOCATE (fsal_wat)
!!$DEALLOCATE (sal_sea)
!!$!
!!$DEALLOCATE (ftfrac)
!!$!
!!$DEALLOCATE (fwc_veg)
!!$DEALLOCATE (fh)
!!$DEALLOCATE (fb)
!!$DEALLOCATE (ftauN)
!!$DEALLOCATE (fw_effL)
!!$DEALLOCATE (fw_effH)
!!$!
!!$DEALLOCATE (fNrh_L)
!!$DEALLOCATE (fNrh_H)
!!$DEALLOCATE (fNrv_L)
!!$DEALLOCATE (fNrv_H)
!!$!
!!$DEALLOCATE (fhrmodel)
!!$!
!!$DEALLOCATE (ftth)
!!$DEALLOCATE (fttv)
!!$!
!!$DEALLOCATE (ftau_veg)
!!$DEALLOCATE (ftau_atm)
!!$DEALLOCATE (ftb_au)
!!$DEALLOCATE (ftb_ad)
!!$DEALLOCATE (ftair)
!!$!
!!$DEALLOCATE (fsnowd)
!!$DEALLOCATE (frsnow)
!!$!
!!$DEALLOCATE (tsoil)
!!$DEALLOCATE (z_lsm)
!!$!
!!$DEALLOCATE (ftheta_inc)  ! PSG: deallocation of pixel-based incidence angle
!!$DEALLOCATE(TB_HV)  ! PSG: deallocating TB multi-dimensional temporal varible
!!$SELECT CASE (CFINOUT)
!!$   CASE ('ascii')
!!$   DEALLOCATE (ndate)
!!$   DEALLOCATE (ntime)
!!$   DEALLOCATE (doy)
!!$   DEALLOCATE (nbpt)
!!$END SELECT
!!$!
