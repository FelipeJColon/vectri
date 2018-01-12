SUBROUTINE setdown
!--------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins AM. 2011, ICTP
! tompkins@ictp.it
!
! setdown subroutine to close shop and tidy up...
!---------------------------------------------------------
  USE netcdf
  USE mo_control
  USE mo_climate
  USE mo_constants
  USE mo_vectri
  USE mo_ncdf_tools

  IMPLICIT NONE

  INTEGER :: latvarid, lonvarid
  INTEGER :: ngonoDimID, &
           & nlarvDimID, &
           & ninfvDimID, &
           & ninfhDimID, &
           & nhostDimID

  !---
  ! 1. close down the ascii and netcdf output files
  !---
  CALL check(NF90_CLOSE(ncidout))

  SELECT CASE(iclim_input)
  CASE(1)
  CASE(2)
  CASE(3)  
    CALL check(NF90_CLOSE(ncid_rain))
    CALL check(NF90_CLOSE(ncid_temp))
    IF (l2d) CALL check(NF90_CLOSE(ncid_pop))
  CASE(4)
  CASE(5)
    CLOSE(infile)
  END SELECT

  write(iounit,*) 'input/output closed - now dump restart'

  !---
  ! 2. dump a restart file
  !---

  ! dump restart file:
!  CALL check(NF90_CREATE(path = input//"restart_vectri.nc", cmode=NF90_CLOBBER, ncid = ncidout))
  CALL check(NF90_CREATE(path = input//"restart_vectri.nc", cmode=or(nf90_clobber,nf90_64bit_offset), ncid = ncidout))


! define the dimensions
  CALL check(NF90_DEF_DIM(ncidout, lon_name, nlon, LonDimID))
  CALL check(NF90_DEF_DIM(ncidout, lat_name, nlat, LatDimID))
  CALL check(NF90_DEF_DIM(ncidout, "ngono", ngono+1, ngonoDimID))
  CALL check(NF90_DEF_DIM(ncidout, "nlarv", nlarv+1, nlarvDimID))
  CALL check(NF90_DEF_DIM(ncidout, "ninfv", ninfv+1, ninfvDimID))

  CALL check(NF90_DEF_DIM(ncidout, "ninfh", ninfh+2, ninfhDimID)) ! immunity +2
  CALL check(NF90_DEF_DIM(ncidout, "nhost", nhost,   nhostDimID))

  ! Define the coordinate variables. They will hold the coordinate
  ! information, that is, the latitudes and longitudes. A varid is
 ! returned for each.
  CALL check(nf90_def_var(ncidout, lon_name, NF90_REAL, londimid, lonvarid) )
  CALL check(nf90_def_var(ncidout, lat_name, NF90_REAL, latdimid, latvarid) )

  ! Assign units attributes to coordinate var data. This attaches a
  ! text attribute to each of the coordinate variables, containing the
  ! units.
  CALL check( nf90_put_att(ncidout, latvarid, units, lat_units) )
  CALL check( nf90_put_att(ncidout, lonvarid, units, lon_units) )

  ! define global attributes
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rundate", now))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "vectri", version))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "command", runcommand))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "comments","vectri restart file"))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "refdate", time_units))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "nday",nday))

  !---------------------------------------------
  ! define the output structures
  ! routine also sets up the diagnostics indices
  ! reuse the nc 2 byte arrays for indices
  !---------------------------------------------
  !                             index      short name title       units
  ndiag2d=0 ! dummy variable

  CALL define_ncdf_output(ncvar_larvae,"larv"  , &
  & "larvae array","m^-2", ndiag2d, (/ nlarvDimID,LonDimId,LatDimID  /))

  CALL define_ncdf_output(ncvar_vector,"vect"  , &
  & "vect array","m^-2", ndiag2d, (/ngonoDimID, ninfvDimId, LonDimID,  LatDimId  /))

  CALL define_ncdf_output(ncvar_pr,"host"  , &
  & "host array","m^-2", ndiag2d, (/ninfhDimID, nhostDimId,  LonDimID, LatDimId  /))
!  & "host array","m^-2", ndiag2d, (/ninfhDimID,  LonDimID, LatDimId  /))

  CALL define_ncdf_output(ncvar_waterfrac,"water_frac"  , &
  & "Water coverage of grid cell","fraction", ndiag2d, (/ LonDimID, LatDimId /) )

  ! *** ADD NEW PROGNOSTIC FIELDS HERE TO RESTART

  ! End define mode.
  CALL check(NF90_ENDDEF(ncidout))

  ! Write the coordinate variable data. This will put the latitudes
  ! and longitudes of our data grid into the netCDF file.
  CALL check( nf90_put_var(ncidout, latvarid, lats) )
  CALL check( nf90_put_var(ncidout, lonvarid, lons) )
  CALL check( nf90_put_var(ncidout, ncvar_larvae(1), rlarv ))
  CALL check( nf90_put_var(ncidout, ncvar_vector(1), rvect ))
  CALL check( nf90_put_var(ncidout, ncvar_pr(1), rhost ))
  CALL check( nf90_put_var(ncidout, ncvar_waterfrac(1), rwaterfrac ))
  CALL check( nf90_close(ncidout))

  RETURN
END
