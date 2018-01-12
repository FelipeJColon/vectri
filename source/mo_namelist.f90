MODULE mo_namelist
!-------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins AM. 2011, ICTP
! tompkins@ictp.it
!
! namelists for input
!---------------------------------------------------------
  USE mo_constants
  USE mo_control
  USE mo_climate
  USE mo_vectri

  IMPLICIT NONE

! program control namelist
  NAMELIST /control/iclim_input,output_file,rundir,ncout_file,start_date,nday,init_file, &
& now,version,runcommand

! climate data namelist 
  NAMELIST /climate/rconst_rain,rconst_temp,rainfile,tempfile

! station data namelist 
  NAMELIST /station/station_file,rpopdensity2010

! region namelist 
  NAMELIST /region/lon,nlon,dlon,lat,nlat,dlat

! constants namelist 
  NAMELIST /constants/neggmn,nlayingmax,rlarv_tmin,rlarv_tmax, &
& rlarv_eggtime,rlarv_pupaetime, rlarv_flushmin,  rlarv_flushtau, &
& nlarv_scheme, rmasslarv_stage4, rbiocapacity, rlarvsurv, &
& rwater_tempoffset, rnobednetuse, rbednettreat, &
& rwaterperm_default, rwaterfrac_rate, &
& rwaterfrac_itau, rwaterfrac_infil125, &
& rwaterfrac_evap126, &
& rwaterfrac_infil130, ipud_vers, rwaterfrac_evap, rwaterfrac_shapep2, &
& rwaterfrac_CN, rwaterfrac_min, rwaterfrac_max,  &
& rbiteratio, rzoophilic_tau, rzoophilic_min, rhostclear, &
& rpthost2vect_I, rpthost2vect_R, rptvect2host, rtgono, dgono, rtsporo, dsporo, &
& nsurvival_scheme, rhost_infect_init, rpop_death_rate, &
& rpopdensity_min, rtemperature_offset, rtemperature_trend, rrainfall_factor, &
& rmigration,rbeta_indoor, &
& rimmune_gain_eira,rimmune_loss_tau, &
& loutput_population, &    
& loutput_rain, &          
& loutput_t2m, &           
& loutput_waterfrac, &     
& loutput_vector, &        
& loutput_larvae, &        
& loutput_lbiomass, &
& loutput_pr, & 
& loutput_prd, &           
& loutput_hbr, &           
& loutput_cspr, &          
& loutput_eir, &           
& loutput_cases, &           
& loutput_immunity, &           
& loutput_vecthostratio, &
& lverbose, &
& nyearspinup, nlenspinup, &
& rbitehighrisk, rvecsurv
  
  END MODULE mo_namelist

