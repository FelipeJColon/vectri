SUBROUTINE initialize
!--------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins AM. 2011, ICTP
! tompkins@ictp.it
!
! initialize subroutine to initialize arrays
!
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

  INTEGER i, j, ncid, ivarid,  nlatcheck, nloncheck
  REAL, ALLOCATABLE :: zdummy2d(:,:)
  REAL :: zfillvalue ! local fill value 

  ! define lat/lon arrays
  DO i=1,nlon
    lons(i)=lon+REAL(i-1)*dlon
  ENDDO
  DO i=1,nlat
    lats(i)=lat+REAL(i-1)*dlat
  ENDDO
  ! set to missing if outside lat band
  DO i=1,nlat 
    IF (ABS(lats(i))>rlatmax) rpopdensity(:,i)=rfillValue
  ENDDO

  ! if needed for stochastic term
  CALL RANDOM_SEED()

  ! initialization options 

  SELECT CASE (TRIM(init_file(1:4)))

  ! 
  ! restart file, exact initialization
  !
  CASE('rest') 

    WRITE(iounit,*) 'opening RESTART condition file '//input//init_file
    CALL check(NF90_OPEN(path=input//init_file,mode=nf90_nowrite,ncid=ncid))

    CALL check(NF90_INQ_DIMID(ncid, "latitude", LatDimID))
    CALL check(NF90_INQUIRE_DIMENSION(Ncid, LatDimID, len = nlatcheck))

    CALL check(NF90_INQ_DIMID(ncid, "longitude", LonDimID))
    CALL check(NF90_INQUIRE_DIMENSION(Ncid, LonDimID, len = nloncheck))
    IF (nlon.ne.nloncheck .or. nlat.ne.nlatcheck) THEN
      WRITE(iounit,*) '*** data dimensions input error *** ',nlon,nloncheck,nlat,nlatcheck,'init_file'
      STOP '*** data dimensions input error *** '
    ENDIF

    !
    ! a. read data for PR initialization
    !
    CALL check(NF90_INQ_VARID(ncid, "host", iVarId))
    CALL check(NF90_GET_VAR(ncid, iVarId, rhost))
    CALL check(NF90_GET_ATT(ncid, iVarId, "_FillValue", zFillValue))
    WHERE(rhost(0,1,:,:)==zfillvalue) rpopdensity=rfillValue

    !
    ! b. read data for larvae initialization
    !
    CALL check(NF90_INQ_VARID(ncid, "larv", iVarId))
    CALL check(NF90_GET_VAR(ncid, iVarId, rlarv))
    !
    ! c. read data for vector initialization
    !
    CALL check(NF90_INQ_VARID(ncid, "vect", iVarId))
    CALL check(NF90_GET_VAR(ncid, iVarId, rvect))
    !
    ! d. read data for water fraction
    !
    CALL check(NF90_INQ_VARID(ncid, "water_frac", iVarId))
    CALL check(NF90_GET_VAR(ncid, iVarId, rwaterfrac))

    ! close the initialization file
    CALL check(NF90_CLOSE(ncid))
    WRITE(iounit,*)'%I: initial conditions read ok'

  CASE ('clim') ! approximate restart from output file...
    !--------------------------------------------------------
    ! initial conditions for PR, vectors, larva and hydrology
    !--------------------------------------------------------
    WRITE(iounit,*) '%I: opening initial condition file '//input//'init_file.nc'
    CALL check(NF90_OPEN(path=input//'init_file.nc',mode=nf90_nowrite,ncid=ncid))
    ! dimensions to make sure input file is correct size.
    CALL check(NF90_INQ_DIMID(ncid, "lat", LatDimID))
    CALL check(NF90_INQUIRE_DIMENSION(Ncid, LatDimID, len = nlatcheck))
    CALL check(NF90_INQ_DIMID(ncid, "lon", LonDimID))
    CALL check(NF90_INQUIRE_DIMENSION(Ncid, LonDimID, len = nloncheck))
    IF (nlon.ne.nloncheck .or. nlat.ne.nlatcheck) THEN
      WRITE(iounit,*) '*** data dimensions input error *** ',nlon,nloncheck,nlat,nlatcheck,'init_file'
      STOP '*** data dimensions input error *** '
    ENDIF

    !
    ! a. read data for PR initialization
    !
    CALL check(NF90_INQ_VARID(ncid, "PR", iVarId))
    CALL check(NF90_GET_VAR(ncid, iVarId, rhost(ninfh,1,:,:)))
    rhost(ninfh,1,:,:)=MAX(rhost(ninfh,1,:,:),0.0)
    rhost(0,1,:,:)=1.0-rhost(ninfh,1,:,:)

    !
    ! b. read data for larvae initialization
    !
    CALL check(NF90_INQ_VARID(ncid, "larvae", iVarId))
    CALL check(NF90_GET_VAR(ncid, iVarId, rlarv(0,:,:)))
    rlarv(0,:,:)=MAX(0.0,rlarv(0,:,:)/FLOAT(nlarv+1))
    DO i=1,nlarv
      rlarv(i,:,:)=rlarv(0,:,:)
    ENDDO

    !
    ! c. read data for vector initialization
    !
    CALL check(NF90_INQ_VARID(ncid, "vector", iVarId))
    CALL check(NF90_GET_VAR(ncid, iVarId, rvect(0,0,:,:)))

    rvect=MAX(0.0,rvect/FLOAT(ngono+1))
    DO j=0,ninfv
      DO i=0,ngono
        rvect(i,j,:,:)=rvect(0,0,:,:)
      ENDDO
    ENDDO

    ALLOCATE(zdummy2d(nlon,nlat))
    CALL check(NF90_INQ_VARID(ncid, "vectinfect", iVarId))
    CALL check(NF90_GET_VAR(ncid, iVarId, zdummy2d(:,:)))

    DO i=0,ngono
      rvect(i,0,:,:)=rvect(i,0,:,:)*(1.0-zdummy2d(:,:))
    ENDDO
    DO j=1,ninfv
      DO i=0,ngono
        rvect(i,j,:,:)=rvect(i,j,:,:)*zdummy2d(:,:)/FLOAT(ninfv)
      ENDDO
    ENDDO
    DEALLOCATE(zdummy2d)

    !
    ! d. read data for puddle
    !
    CALL check(NF90_INQ_VARID(ncid, "water_frac", iVarId))
    CALL check(NF90_GET_VAR(ncid, iVarId, rwaterfrac(:,:)))
    rwaterfrac=MAX(0.0,MIN(rwaterfrac,1.0))

    ! close the initialization file
    CALL check(NF90_CLOSE(ncid))

  CASE DEFAULT  ! idealized fixed initial conditions 
    WRITE(iounit,*)'modelling initializing from artificial conditions'

    rzoophilic=0.0
    rbitezoo=0.0

    rvect=0.0
    rlarv=0.0
    rhost=0.0
    rpuddle=0.0 
    rwaterfrac=rwaterfrac_min

    ! initialization for 
    ! reservoir of malaria for initial conditions...
    ! can be overwritten later
    rhost(0,1,:,:)=1.0-rhost_infect_init ! set % that start with parasite
    rhost(ninfh,1,:,:)=rhost_infect_init
    rvect(:,0,:,:)=100*rhost_infect_init*rvect_min          ! fixed density of vectors
    rvect(:,ninfv,:,:)=10*rhost_infect_init*rvect_min          ! fixed density of vectors

  END SELECT


  RETURN
END SUBROUTINE initialize
