SUBROUTINE setup
!--------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins AM. 2011, ICTP
! tompkins@ictp.it
!
! setup subroutine to read in forcing data and setup output files
!
! *** THIS ROUTINE IS A MESS AND TOO LONG ***
!     WILL BE SUBDIVIDED INTO SUBTASKS SOON
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
  INTEGER :: ix,iy,i

  !WRITE(iounit,*) 'reading namelists'
  CALL read_namelists

  !
  ! allocate arrays
  !
  ALLOCATE(lons(nlon))
  ALLOCATE(lats(nlat))
  ALLOCATE(rwaterfrac(nlon,nlat))   ! diagnostic 
  ALLOCATE(rpuddle(nlon,nlat)) ! prognostic
  ALLOCATE(rwaterperm(nlon,nlat))   ! diagnostic
  ALLOCATE(rzoophilic(nlon,nlat))
  ALLOCATE(rbitezoo(nlon,nlat))
  ALLOCATE(rvect(0:ngono,0:ninfv,nlon,nlat))
  ALLOCATE(rlarv(0:nlarv,nlon,nlat))   

  ! v1.4.0 add one to ninfh for immunity:
  ALLOCATE(rhost(0:ninfh+1,nhost,nlon,nlat))  

  ALLOCATE(rtemp(nlon,nlat))
  ALLOCATE(rrain(nlon,nlat))
  ALLOCATE(rpopdensity(nlon,nlat))
  ALLOCATE(rwateroccupancy(nlon,nlat))


  ! safety check on the spin up:  
  IF (nyearspinup>0 .AND. nday<nlenspinup) THEN 
     WRITE(iounit,*) 'spinup period reset to integration length'
     nlenspinup=nday
  ENDIF

  !WRITE(iounit,*) 'open input'
  CALL open_input

  !--------------------------------------
  ! initialize prognostic data constructs
  !--------------------------------------
  !WRITE(iounit,*) 'initialize model'
  CALL initialize

  !------------------------------------
  ! migration scaled to be per timestep
  !------------------------------------
  ! rmigration=rmigration*dt/rdaysinyear

  !   Initial values for vectors/larvae...
  DO ix=1,nlon
  DO iy=1,nlat
    IF (rpopdensity(ix,iy)<rpopdensity_min) THEN  ! missing population data means a sea point
      rhost(:,:,ix,iy)=rfillvalue ! default value to prevent divide by zero error
      rvect(:,:,ix,iy)=rfillvalue ! default value to prevent divide by zero error
      rlarv(:,ix,iy)=rfillvalue ! default value to prevent divide by zero error
      rwaterfrac(ix,iy)=rfillvalue ! default value to prevent divide by zero error
      rpopdensity(ix,iy)=rfillvalue
    ENDIF
  ENDDO 
  ENDDO

  ! mass relationship of larvae - we assume a stage 4 size as Bomblies and 
  ! we will assume a linear mass increase with age unless find other reference 
  ! this linear growth rate is very close to that assumed by Bomblies.   
  DO i=0,nlarv
    rmasslarv(i)=FLOAT(i)*rmasslarv_stage4/FLOAT(nlarv)
  ENDDO

  ! define the permanent water fraction 
  ! IF (lperm_sentinel) THEN
  !   read in the sentinel map...
  ! ELSE 
      rwaterperm(:,:)=rwaterperm_default
  ! ENDIF

  !
  ! tranmission probabilities as a function of bite number (for speed)
  !
  rpdfvect2host=0.0
  DO i=1,nbitepdf
    rpdfvect2host(i)=1.0-(1.0-rptvect2host)**REAL(i)
  ENDDO

  ! initialize the immunity gain function
  ! conversion of immune rate into tau parameter
  ! assumes 95% are immune at this value
  rimmune_gain_tau=-rimmune_gain_eira/LOG(0.05)

  ! convert annual death rate to daily 
  rpop_death_rate_daily=rpop_death_rate/365.0

  ! set the rwaterfrac_S value for puddle 130
  rwaterfrac_S=25400/rwaterfrac_CN-254

  ! place the larval progression rates into an array
  rlarvmature(1,:)=(/0.0,rlarv_ermert/)  ! temperature independent
  rlarvmature(2,:)=(/rlarv_jepson,-rlarv_tmin*rlarv_jepson/)
  rlarvmature(3,:)=(/rlarv_bayoh,-rlarv_tmin*rlarv_bayoh/)
  rlarvmature(4,:)=(/rlarv_craig,-rlarv_tmin*rlarv_craig/)

  ! -------------
  ! safety checks
  ! -------------
  nxdg=MAX(MIN(nxdg,nlon),1)
  nydg=MAX(MIN(nydg,nlat),1)

  !write(iounit,*),'opening output'
  CALL open_output

CONTAINS
!-----------------------------------------------
  SUBROUTINE read_namelists

!-------------------------------
! 1. NAMELISTS
!-------------------------------
  !write(iounit,*) 'file ',input//'vectri.namelist'
  OPEN(8,file=input//'vectri.namelist',status='OLD')
  READ(8,nml=control)
  READ(8,nml=climate)
  READ(8,nml=station)
  READ(8,nml=region)
  READ(8,nml=constants)
  CLOSE(8) ! namelist file

!-------------------------------
! 2. OUTPUT NAMELIST DATA
!-------------------------------

! ------
! output
! ------
  IF (output_file=='screen') THEN
    iounit=6
  ELSE
    iounit=7
    OPEN(iounit,FILE=output//'vectri.out')
  ENDIF

  WRITE(iounit,*) '-------------- VECTRI '//version//' ---------------'
  WRITE(iounit,*) 'run date ',now
  WRITE(iounit,*)

  WRITE(iounit,*)' integration days: ',nday
!  WRITE(iounit,*)' timestep: ',dt
!  WRITE(iounit,*)' integration steps: ',nstep
  WRITE(iounit,*) 
  WRITE(iounit,*)' advection method: ',nnumeric
  WRITE(iounit,*) 


  WRITE(iounit,*)' Climate namelist'
  WRITE(iounit,*)' climate method (constant/point/gridded): ',iclim_input
  SELECT CASE(iclim_input)
  CASE(1)
    WRITE(iounit,*)' constant input selected'
    WRITE(iounit,*)' constant rain(mm/day): ',rconst_rain 
    WRITE(iounit,*)' constant temperature(C): ',rconst_temp 
  CASE(2)
!stop 'write here for netcdf'
    WRITE(iounit,*)' GTS netcdf station data input selected '
    WRITE(iounit,*)' station file is ',station_file
  CASE(3)
    WRITE(iounit,*)' gridded data selected '
    WRITE(iounit,*)' data source ',station_file
    WRITE(iounit,*)' base year 2010 population density (/km2) ',rpopdensity2010
  CASE(4)
    WRITE(iounit,*)' station data input selected '
    WRITE(iounit,*)' station file is ',station_file
  CASE(5)
    WRITE(iounit,*)' ADDIS ASCII station data input selected '
  CASE DEFAULT
    WRITE(iounit,*)' Please select data input type from 1 to 5'
    STOP
  END SELECT
  
  RETURN
END SUBROUTINE read_namelists

END SUBROUTINE setup


