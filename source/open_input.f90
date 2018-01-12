SUBROUTINE open_input
!--------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins AM. Aug 2012, ICTP
! tompkins@ictp.it
!
! open and read the input files... added to V1.3
!
!---------------------------------------------------------
  USE netcdf
  USE mo_control
  USE mo_climate
  USE mo_constants
  USE mo_vectri
  USE mo_namelist
  USE mo_ncdf_tools

  IMPLICIT NONE

  ! local variables
  INTEGER :: ix,iy,i, istatus=0
  CHARACTER(LEN=10) :: dummy

  ! netcdf vars - only needed locally
  INTEGER :: iVarId 
  INTEGER :: ndaycheck
  REAL :: tlats(nlat), tlons(nlon) ! temporary lat and lon arrays

  ! logical to state if there are 2D arrays.
  l2d=((nlon>1.or.nlat>1).and.iclim_input==3)

!---------------------------------------------------------------------------
!  meteorology variables read in from analysis, obs or forecast, or constant
!---------------------------------------------------------------------------
  SELECT CASE(iclim_input)

  !-----------------
  ! constant forcing
  !-----------------
  CASE(1)
    nstep=nday/NINT(dt) ! length of run in timesteps
!    ALLOCATE(rrh(nlon,nlat,nday))
    ALLOCATE(ndate(nday)) 
    nrun=(nday+nyearspinup*INT(nlenspinup))/NINT(dt)

    ndaydiag=360

    DO ix=1,nlon
    DO iy=1,nlat
      DO i=1,nday
        ndate(i)=i
      ENDDO
!     rtemp(ix,iy,:)=rconst_temp
!     rrain(ix,iy,:)=rconst_rain
    ENDDO
    ENDDO

    ! population density fixed
    rpopdensity(:,:)=rpopdensity2010

    ! 1.e-6 is to convert popdensity into SI units (km->m) !!!
    rpopdensity(:,:)=rpopdensity(:,:)*1.e-6

  CASE(2)

    STOP 'code not ready for generalized single point netcdf yet'
    
  !------------------------------------------------------------------------
  ! gridded input - assumes the input files have the correct timestamps
  !                 i.e. just reads the file blindly - the lat/lon info is 
  !                 important though...
  !------------------------------------------------------------------------
  CASE(3)
    !-----------------------------------------------------------------------------
    ! rainfall
    ! lat and lon now read from rainfall file and checked against temperature file
    !-----------------------------------------------------------------------------
    WRITE(iounit,*)'%I: opening rainfall ',rainfile
    CALL check(NF90_OPEN(path=rainfile,mode=nf90_nowrite,ncid=ncid_rain))

    ! for 2D field check dimensions to make sure input file is correct size.
    ! NOTE! lats and lons now read in from rainfile...
    IF (l2d) CALL check_latlon(ncid_rain,tlats,tlons)

    ! read in number of timesteps
    CALL check(NF90_INQ_DIMID(ncid_rain, "time", timeDimID))

    ! read data for rainfall, check standard names
    DO i=1,SIZE(rainvar)
      istatus=NF90_INQ_VARID(ncid_rain, TRIM(rainvar(i)), VarId_rain)
      IF(istatus==nf90_noerr)EXIT
    ENDDO

    IF (istatus /= nf90_noerr) THEN
      WRITE(IOUNIT,*)'%WARNING: rain name not found, using brute force approach'
      VarId_rain=2 ! not a predetermined name so read first variable by brute force
    ENDIF

    ! if attribute not there then default to fillvalue
    istatus=NF90_GET_ATT(ncid_rain, VarId_rain, "_FillValue", rrain_FillValue)
    IF (istatus /= nf90_noerr) rrain_FillValue=rfillvalue
    IF (lclim_chunk) CALL check(NF90_GET_VAR(ncid_rain, VarId_rain, rrain ))

    !------------
    ! temperature
    !------------
    WRITE(iounit,*)'%I: opening temperature ',tempfile
    CALL check(NF90_OPEN(path=tempfile,mode=nf90_nowrite,ncid=ncid_temp))
    ! dimensions to make sure input file is correct size.
    IF (l2d) call check_latlon(ncid_temp,tlats,tlons)
    ! Use simple checksum to make sure temp and rain files match
!    IF (SUM(lats)/=SUM(tlats).or.SUM(lons)/=SUM(tlons)) &
!    &  STOP 'ERROR lat/lon arrays different between climate files'
    
    ! data for temperature, try to open with standard name first
    ! read data for rainfall, check standard names
    DO i=1,SIZE(tempvar)
      istatus=NF90_INQ_VARID(ncid_temp, TRIM(tempvar(i)), VarId_temp)
      IF(istatus==nf90_noerr)EXIT
    ENDDO
    IF (istatus/=nf90_noerr) THEN
      WRITE(IOUNIT,*)'+++ WARNING: temperature name not found, using brute force approach'
      VarId_temp=2 ! brute force approach
    ENDIF

    IF (lclim_chunk) CALL check(NF90_GET_VAR(ncid_temp, VarId_temp, rtemp ))

    ! if attribute not there then default to fillvalue
    istatus=NF90_GET_ATT(ncid_temp, VarId_temp, "_FillValue", rtemp_FillValue)
    IF (istatus /= nf90_noerr) rtemp_FillValue=rfillvalue

    ! allocate the fields for dates (has to be after above, since nday read from rainfile
    ALLOCATE(ndate(nday)) 
    nrun=(nday+nyearspinup*INT(nlenspinup))/NINT(dt) ! NEED TO CHANGE - 3D=nday

    ! dates
    ndaydiag=1
    DO i=1,nday    
      ndate(i)=i-1
    ENDDO
    nstep=nday/NINT(dt) ! length of run in timesteps

    !----------------
    ! population file
    !----------------
    WRITE(iounit,*)'%I: opening population ',ncpopufile
    CALL check(NF90_OPEN(path=ncpopufile,mode=nf90_nowrite,ncid=ncid_pop))
    ! dimensions to make sure input file is correct size.
    CALL check_latlon(ncid_pop,tlats,tlons)

    !---------------------
    ! read population data
    !---------------------
    CALL check(NF90_INQ_VARID(ncid_pop, "population", iVarId))
    CALL check(NF90_GET_ATT(ncid_pop, iVarId, "_FillValue", rpopdensity_FillValue))

    ! simple fixed population with fixed growth - otherwise will link to population MODEL! 
    IF (kpop_option==1) CALL check(NF90_GET_VAR(ncid_pop, iVarId, rpopdensity(:,:)))

    !----------------
    ! LUC occupancy file
    !----------------
    IF (.false.) THEN
      WRITE(iounit,*)'%I: opening land use occupancy ',nclucfile
      CALL check(NF90_OPEN(path=nclucfile,mode=nf90_nowrite,ncid=ncid_luc))
      ! dimensions to make sure input file is correct size.
      CALL check_latlon(ncid_pop,tlats,tlons)

    !---------------------
    ! read data
    !---------------------
      CALL check(NF90_INQ_VARID(ncid_luc, "Band1", iVarId))
      ! CALL check(NF90_GET_ATT(ncid_pop, iVarId, "_FillValue", rpopdensity_FillValue))

    ! simple fixed population with fixed growth - otherwise will link to population MODEL! 
      CALL check(NF90_GET_VAR(ncid_luc, iVarId, rwateroccupancy(:,:)))
    ELSE  ! single grid point OR NEED TO ADD SWITCH FOR CONSTANT VALUE
      rwateroccupancy=1.0 ! or fixed value
    ENDIF !l2d


    ! -----------------------------
    ! scale factors and processing!
    ! -----------------------------
    WHERE(rpopdensity<0.0)rpopdensity=rfillValue
    
    ! 1.e-6 is to convert popdensity into SI units (km->m) !!!
    WHERE(rpopdensity/=rfillValue)rpopdensity=rpopdensity*1.e-6

    ! all sparsely populated point set to missing.
    WHERE(rpopdensity<rpopdensity_min)rpopdensity=rfillValue

    ! all negative water scal fac set to zero
    WHERE(rwateroccupancy<0.0)rwateroccupancy=0.0

  !-------------------------------
  ! point station data from Volker
  !-------------------------------
  CASE(4)
    ndaydiag=360
    ALLOCATE(ndate(nday)) 
    nrun=(nday+nyearspinup*INT(nlenspinup))/NINT(dt) ! NEED TO CHANGE - 3D=nday
    OPEN(infile,file=station_file,action='read')      
    READ(infile,*) dummy
    ! population density fixed
    ! 1.e-6 is to convert popdensity into SI units (km->m) !!!
    rpopdensity(:,:)=rpopdensity2010*1.e-6

  !------------------------------------
  ! ASCII point station data from Gizaw 
  !------------------------------------
  CASE(5)
    ndaydiag=365
    ALLOCATE(ndate(nday)) 
    nrun=(nday+nyearspinup*INT(nlenspinup))/NINT(dt) ! NEED TO CHANGE - 3D=nday
    OPEN(infile,file=input//station_file,action='read')      
    ! population density fixed
    ! 1.e-6 is to convert popdensity into SI units (km->m) !!!
    rpopdensity=rpopdensity2010*1.e-6
  END SELECT

  DO i=1,nday    
    ndate(i)=i-1
  ENDDO

CONTAINS
!---------------------------------------------------------

SUBROUTINE check_latlon(ncid,tlats,tlons)

  USE mo_control
  USE mo_vectri

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ncid
  REAL, INTENT(INOUT) :: tlats(nlat),tlons(nlon)
  INTEGER :: nlatcheck,nloncheck,istatus

  istatus=NF90_INQ_DIMID(ncid, "lat", LatDimID)
  IF (istatus /= nf90_noerr) CALL check(NF90_INQ_DIMID(ncid, "latitude", LatDimID))
  CALL check(NF90_INQUIRE_DIMENSION(ncid, LatDimID, len = nlatcheck))

  istatus=NF90_INQ_DIMID(ncid, "lon", LonDimID)
  IF (istatus /= nf90_noerr) CALL check(NF90_INQ_DIMID(ncid, "longitude", LonDimID))
  CALL check(NF90_INQUIRE_DIMENSION(ncid, LonDimID, len = nloncheck))

! CALL check(NF90_INQ_DIMID(ncid, "time", timeDimID))
! CALL check(NF90_INQUIRE_DIMENSION(ncid, timeDimID, len = ndaycheck))

! here we check that the temperature file has the same number of days as the rain file
  IF (nlon.ne.nloncheck .or. &
    &   nlat.ne.nlatcheck ) THEN !.or. &
!    &   nday.ne.ndaycheck) THEN
    WRITE(iounit,*) '*** data dimensions input error *** ',nlon,nloncheck,nlat,nlatcheck,nday,ndaycheck, tempfile
    STOP '*** data dimensions input error *** '
  ENDIF

! read in the lat and lon arrays to pass back.
  CALL check(NF90_GET_VAR(ncid, LatDimID, tlats ))
  CALL check(NF90_GET_VAR(ncid, LonDimID, tlons ))

RETURN
END SUBROUTINE check_latlon

END SUBROUTINE open_input
