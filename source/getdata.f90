SUBROUTINE getdata(iday)
!--------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins AM. 2012, ICTP
! tompkins@ictp.it
!
! getdata subroutine to read in climate and population slice
!
!
!---------------------------------------------------------
  USE netcdf
  USE mo_control
  USE mo_climate
  USE mo_constants
  USE mo_vectri
  USE mo_ncdf_tools

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: iday

  REAL    :: ztempmax,ztempmin
  INTEGER :: iyyyy,imm,idd,idate

  !---------------------------------------------------------------------------
  !  meteorology variables read in from analysis, obs or forecast, or constant
  !---------------------------------------------------------------------------
  SELECT CASE(iclim_input)

  !-----------------
  ! constant forcing
  !-----------------
  CASE(1)
    rtemp=rconst_temp
    rrain=rconst_rain

  ! NETCDF station data
  CASE(2)
    STOP

  !------------------------------------------------------------------------
  ! gridded input - assumes the input files have the correct timestamps
  !                 i.e. just reads the file blindly - the lat/lon info is 
  !                 important though...
  !------------------------------------------------------------------------
  CASE(3)
    IF (.not. lclim_chunk) THEN
      CALL check(NF90_GET_VAR(ncid_rain, VarId_rain, rrain, start=(/1,1,iday/),count=(/nlon,nlat,1/)))
      CALL check(NF90_GET_VAR(ncid_temp, VarId_temp, rtemp, start=(/1,1,iday/),count=(/nlon,nlat,1/)))
    ENDIF

    WHERE(rtemp==rtemp_Fillvalue) rtemp=rfillvalue ! set to negative default

    ! scalfac for TRMM required
    IF (l2d) THEN
      WHERE(rrain==rrain_Fillvalue) 
        rrain=rfillvalue
      ENDWHERE
    ENDIF
  CASE(4)
    READ(infile,*) idate,rconst_temp,rconst_rain 
    ndate(iday)=idate
    rtemp=rconst_temp
    rrain=rconst_rain
  CASE(5)
    IF (iday==1) REWIND(infile)
    READ(infile,*) iyyyy,imm,idd,rconst_rain,ztempmin,ztempmax 
    rconst_temp=(ztempmin+ztempmax)/2.0
    ! quality control
    IF (rconst_temp>rqualtemp_max) rconst_temp=rqualtemp_max
    IF (rconst_temp<rqualtemp_min) rconst_temp=rqualtemp_min

    ndate(iday)=iday !idate
    rtemp=rconst_temp
    rrain=rconst_rain
  END SELECT

  ! ----------------------------------- 
  ! Unit conversion and quality control
  ! ----------------------------------- 
  ! ERAI convert K to C if T_mean >240 (either warm or in K!) 
  WHERE(rtemp>240.) rtemp=rtemp-r0CtoK 

  ! quality control for bad station data
  WHERE(rtemp>60.) rtemp=20.0
  WHERE(rtemp<-10.) rtemp=20.0

END SUBROUTINE getdata

