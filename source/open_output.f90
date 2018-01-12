SUBROUTINE open_output
!--------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins AM. aug 2012, ICTP
! tompkins@ictp.it
!
! opens the output files - added to v1.3
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

  ! netcdf vars - only needed locally
  INTEGER :: latvarid, lonvarid, timevarid

  ! output file details
  CHARACTER (len = 8)  :: str_date
!----------------------------
! 3. OUTPUT FILES
!----------------------------
  !---------------------------
  ! OUTPUT NETCDF FILE
  !---------------------------
  WRITE(str_date,'(I8)') start_date
  time_units="days since "//str_date(1:4)//"-"//str_date(5:6)//"-"//str_date(7:8)//" 00:00 UTC"

  CALL check(NF90_CREATE(path = ncout_file, cmode=NF90_CLOBBER, ncid = ncidout))

! define the dimensions
  CALL check(NF90_DEF_DIM(ncidout, lon_name, nlon, LonDimID))
  CALL check(NF90_DEF_DIM(ncidout, lat_name, nlat, LatDimID))
  CALL check(NF90_DEF_DIM(ncidout, time_name, nf90_unlimited, timeDimID))

  ! Define the coordinate variables. They will hold the coordinate
  ! information, that is, the latitudes and longitudes. A varid is
  ! returned for each.
  CALL check(nf90_def_var(ncidout, lon_name, NF90_REAL, londimid, lonvarid) )
  CALL check(nf90_def_var(ncidout, lat_name, NF90_REAL, latdimid, latvarid) )
  CALL check(nf90_def_var(ncidout, time_name, NF90_REAL, timedimid, timevarid) )

  ! Assign units attributes to coordinate var data. This attaches a
  ! text attribute to each of the coordinate variables, containing the
  ! units.
  CALL check( nf90_put_att(ncidout, latvarid, units, lat_units) )
  CALL check( nf90_put_att(ncidout, lonvarid, units, lon_units) )
  CALL check( nf90_put_att(ncidout, timevarid, units, time_units) )

  ! define global attributes
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "origin","Felipe J Colon-Gonzalez, University of East Anglia, Norwich, UK"))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "web","www.uea.ac.uk/environmental-sciences/people/profile/f-colon"))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "email","f.colon-AT-uea.ac.uk"))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "license","free for research educational purposes, no third party distribution"))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "date", now))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "version", version))


  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rainfile", TRIM(rainfile)))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "tempfile", TRIM(tempfile)))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "neggmn", neggmn))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "nlayingmax",nlayingmax))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rlarv_tmin",rlarv_tmin))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rlarv_tmax",rlarv_tmax))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rlarv_eggtime",rlarv_eggtime))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rlarv_pupaetime",rlarv_pupaetime))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rlarv_flushmin",rlarv_flushmin))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rlarv_flushtau",rlarv_flushtau))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rlarv_ermert",rlarv_ermert))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rlarv_jepson",rlarv_jepson))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rlarv_bayoh",rlarv_bayoh))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rbeta_indoor",rbeta_indoor))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rmasslarv_stage4",rmasslarv_stage4))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rbiocapacity",rbiocapacity))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rlarvsurv",rlarvsurv))

  ! output depends on hydrology version
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rwaterperm_default", rwaterperm_default))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rwaterfrac_rate", rwaterfrac_rate))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "ipud_vers", ipud_vers))

  SELECT CASE(ipud_vers)
  CASE(125)
    CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rwaterfrac_itau", rwaterfrac_itau))
    CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rwaterfrac_infil125", rwaterfrac_infil125))
  CASE(126)
    CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rwaterfrac_evap126", rwaterfrac_evap126))
    CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rwaterfrac_max", rwaterfrac_max))
  CASE(130)
    CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rwaterfrac_CN", rwaterfrac_CN))
    CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rwaterfrac_min", rwaterfrac_min))
    CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rwaterfrac_max", rwaterfrac_max))
    CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rwaterfrac_evap", rwaterfrac_evap))
    CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rwaterfrac_infil130", rwaterfrac_infil130))
    CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rwaterfrac_shapep2", rwaterfrac_shapep2))
  END SELECT

  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rwater_tempoffset",rwater_tempoffset ))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rnobednetuse",rnobednetuse))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rbednettreat",rbednettreat))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rbiteratio",rbiteratio))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rbitehighrisk",rbitehighrisk))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rzoophilic_tau",rzoophilic_tau))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rzoophilic_min",rzoophilic_min))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rhostclear",rhostclear))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rpthost2vect_I",rpthost2vect_I))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rpthost2vect_R",rpthost2vect_R))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rptvect2host",rptvect2host))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rtgono",rtgono))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "dgono",dgono))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rtsporo",rtsporo))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "dsporo",dsporo))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "nlarv_scheme",nlarv_scheme))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rvecsurv",rvecsurv))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "nsurvival_scheme",nsurvival_scheme))
!  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rtlarsurvmax",rtlarsurvmax))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rmar1",rmar1))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rmar2",rmar2))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rtvecsurvmin",rtvecsurvmin))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rtvecsurvmax",rtvecsurvmax))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rhost_infectd",rhost_infectd))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rhost_detectd",rhost_detectd))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rimmune_gain_eira",rimmune_gain_eira))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rimmune_loss_tau",rimmune_loss_tau))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rpop_death_rate",rpop_death_rate))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rpopdensity_min",rpopdensity_min))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "nyearspinup",nyearspinup))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "nlenspinup",nlenspinup))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "dt",dt))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "nnumeric",nnumeric))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rtemperature_offset",rtemperature_offset))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rtemperature_trend",rtemperature_trend))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rrainfall_factor",rrainfall_factor))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rmigration",rmigration))
  CALL check( NF90_PUT_ATT(ncidout, NF90_GLOBAL, "rhost_infect_init",rhost_infect_init))

  !---------------------------------------------
  ! *** define the 3d structures
  ! routine also sets up the diagnostics indices
  !---------------------------------------------
  !                             index      short name title       units
  ndiag2d=0

  IF(loutput_population) &
  & CALL define_ncdf_output(ncvar_popdensity,"population", &
  & "population density","m^-2",  ndiag2d, (/LonDimId, LatDimID, timeDimID/))

  IF(loutput_rain) &
  & CALL define_ncdf_output(ncvar_rain,"rain"    ,"Rainfall","mm/day", ndiag2d, (/ LonDimId, LatDimID, timeDimID /) )

  IF(loutput_t2m) &
  & CALL define_ncdf_output(ncvar_t2m,"t2m"    ,"Temperature","degrees C", ndiag2d, (/ LonDimId, LatDimID, timeDimID /) )

  IF(loutput_waterfrac) &
  & CALL define_ncdf_output(ncvar_waterfrac,"waterfrac" , &
  & "Water coverage of grid cell","fraction", ndiag2d, (/ LonDimId, LatDimID, timeDimID /) )

  IF(loutput_vector) &
  & CALL define_ncdf_output(ncvar_vector,"vector"  ,"Vector Number","density m^-2", ndiag2d,(/LonDimId, LatDimID, timeDimID/))

  IF(loutput_larvae) &
  & CALL define_ncdf_output(ncvar_larvae,"larvae"  ,"larvae Number","density m^-2", ndiag2d, (/LonDimId, LatDimID, timeDimID/))

  IF(loutput_lbiomass) &
  & CALL define_ncdf_output(ncvar_lbiomass,"lbiomass"  , &
  & "larvae biomass","mg m^-2",  ndiag2d, (/LonDimId, LatDimID, timeDimID/))

  IF(loutput_pr) &
  & CALL define_ncdf_output(ncvar_pr,"PR"  , &
  & "Parasite Ratio - Proportion of hosts infected","frac",  ndiag2d, (/ LonDimId, LatDimID, timeDimID /) )
!  CALL define_ncdf_output(ncvar_vectinfect,"vectinfect"  , &
!   & "Proportion of vectors infected","frac",  ndiag2d, (/ LonDimId, LatDimID, timeDimID /) )

  IF(loutput_prd) &
  & CALL define_ncdf_output(ncvar_prd,"PRd", &
  &  "detectable parasite ratio - (after day 10)", &
  & "frac",  ndiag2d, (/ LonDimId, LatDimID, timeDimID /) )

  IF(loutput_hbr) &
  & CALL define_ncdf_output(ncvar_hbr,"hbr"  ,"HBR - Human bite rate", &
  & "bites per day per person",  ndiag2d, (/ LonDimId, LatDimID, timeDimID /) )

  IF(loutput_cspr) & 
  & CALL define_ncdf_output(ncvar_cspr,"cspr" , &
  &  "circumsporozoite protein rate - proportion of infective vector","frac",  ndiag2d, (/ LonDimId, LatDimID, timeDimID /))

  IF(loutput_eir) & 
  & CALL define_ncdf_output(ncvar_eir,"eir" , &
  &  "Entomological Inoculation Rate","infective bites per day per person",  ndiag2d, (/ LonDimId, LatDimID, timeDimID /) )

  IF(loutput_cases) & 
  & CALL define_ncdf_output(ncvar_cases,"cases" , &
  &  "Number of new cases","fraction",  ndiag2d, (/ LonDimId, LatDimID, timeDimID /) )

  IF(loutput_immunity) & 
  & CALL define_ncdf_output(ncvar_immunity,"immunity" , &
  &  "Immune propulation","fraction",  ndiag2d, (/ LonDimId, LatDimID, timeDimID /) )

  IF(loutput_vecthostratio) & 
  & CALL define_ncdf_output(ncvar_vecthostratio,"vecthostratio"  , &
  & "vector to host ratio","frac",  ndiag2d, (/ LonDimId, LatDimID, timeDimID /) )

  ! End define mode.
  CALL check(NF90_ENDDEF(ncidout))

  ! Write the coordinate variable data. This will put the latitudes
  ! and longitudes of our data grid into the netCDF file.
  CALL check( nf90_put_var(ncidout, latvarid, lats) )
  CALL check( nf90_put_var(ncidout, lonvarid, lons) )
  CALL check( nf90_put_var(ncidout, timevarid, ndate) )

  !  IF (loutput_rain) CALL check( nf90_put_var(ncidout, ncvar_rain(1), rrain ))
  !  IF (loutput_t2m) CALL check( nf90_put_var(ncidout, ncvar_temp(1), rtemp ))

  ! NOTE we rescale to non-SI units for NETCDF output for ease of use.
  IF(loutput_population) CALL check( nf90_put_var(ncidout, ncvar_popdensity(1), rpopdensity ))

  !------------------
  ! Set up the arrays
  !------------------
  ! at end of list ndiag contains the number of diagnostics... so we can allocate the diags array
  ! NOTE: ***NEVER*** use rdiag in any calculations of cross variables as the user can switch these
  !       off.  rdiag is for diagnostic storage ONLY! Any other use is unsafe.
  ALLOCATE(rdiag2d(nlon,nlat,ndiag2d))
  rdiag2d=0.0

  RETURN
END SUBROUTINE open_output
