MODULE mo_control
!--------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins AM. 2011, ICTP
! tompkins@ictp.it
!
! run control parameters for model with initial values
!---------------------------------------------------------

  USE mo_constants

  IMPLICIT NONE

!
! these will be user definable soon
!  
!
! run control
!
  INTEGER :: iounit=7
  CHARACTER (len=500) :: output_file, ncout_file, rundir
  CHARACTER(len=*), PARAMETER :: input='./input/'
  CHARACTER(len=*), PARAMETER :: output='./output/'

! these...
  CHARACTER (len = 50) :: now,version
  CHARACTER (len = 500) :: runcommand

! ASCII file output in additional to NETCDF - defunct from v1.3.3
  LOGICAL :: lascii=.false.
  LOGICAL :: lverbose=.false.
!
! Larvae scheme 
!   1: constant as in Ermert 2010
!   2: Jepson 1947 as reported by HM04
!   3: BL2003 (approximated linearly)
!   4: Craig et al. (1999) - close to BL2003
!
  INTEGER :: nlarv_scheme=4

! vector temperature dependent survival scheme in use (emert 2010)
  ! 1= MartinsI
  ! 2= MartinsII
  ! 3= Bayoh
  INTEGER :: nsurvival_scheme=2 

  REAL, PARAMETER :: dt =24./24.       ! timestep in days, ~HAS TO BE DAILY FOR MOMENT =1

! Time constants
  REAL :: rdaysinyear=365.0    ! number of days in year (LMM uses 360)

! model spin up - nyearspinup now refers to the number of spin up "loops" - and the loop length
! is set by nlenspinup.  Examples:
! Repeat the first year 3 times:          nyearspinup=3  nlenspinup=365
! Repeat the first 60 days, 10 times:   nyearspinup=10 nlenspinup=60
  INTEGER  :: nyearspinup=0    ! Number of loops (no longer fixed to a year, name same for compatibility)
  INTEGER  :: nlenspinup=365   ! spinup period of one loop (units of days)

! move derived parameters to setup.f90
  INTEGER :: nstep ! =nday/dt ! length of run in timesteps ASSUME DT=1 for moment
  
! numerics parameters
!
! numeric scheme for advection
! 1:  - box category HM04
! 2:  - split box (inverse SL)
! 3:  - semi langrangian linear interpolation
! 4:  - semi langrangian cubic interpolation
! 5:  - semi langrangian GTD cubic interpolation
! 6:  - semi langrangian cubic spline interpolation
  INTEGER :: nnumeric=2

! array resolution definition constants
  INTEGER, PARAMETER :: nlarv=25  ! resolution of larvae development
! note: NGONO should be an EVEN NUMBER:
  INTEGER, PARAMETER :: ngono=25  ! resolution of mosquito gonotropic development
  INTEGER, PARAMETER :: ninfv=25  ! resolution of vector parasite development 

  INTEGER, PARAMETER :: ninfh=NINT(rhost_infectd/dt)   ! resolution of host parasite development 

  ! Array size of eir pdf distributions
  ! No point setting this parameter to values > 20 since infection probability is already 0.999 
  INTEGER, PARAMETER :: nbitepdf=5  

! trickle minimum vector concentration due to accidental import or advection
  REAL :: rvect_min=1.e-6/ngono  ! 1.3.2 bug correct, remove scale by ninfect

! index of point to dump quick timeseries diagnostics
  INTEGER :: nxdg=427,nydg=226 

! netcdf and input file indices
  INTEGER :: ncid_temp, ncid_rain, ncid_pop, ncid_luc, infile=8    ! file id
  INTEGER :: varid_temp, varid_rain, varid_pop ! var  id
!  LOGICAL :: ltemp_latreverse=.false. ! reverse the lats?
!  LOGICAL :: lrain_latreverse=.false.

! metdata files
  INTEGER  :: nrainsource=0 ! will be read in from file... 0 for model/forecast.
! In revamped v1.4.0 the input/output will need to be cleaned up...
  CHARACTER, PARAMETER :: ncpopufile*100=input//'popufile.nc'
  CHARACTER, PARAMETER :: nclucfile*100=input//'lucfile.nc'
  INTEGER :: kpop_option=1
  CHARACTER*100 :: init_file='none'
  REAL :: rfillvalue=-9999.0 ! missing value

  ! climate control
  ! lclim_chunk controls how climate info read in 
  ! true - all in one go, heavy on memory but much faster
  ! false - saves memory but model runs much slower.
  LOGICAL :: lclim_chunk=.false.  !true is currently bugged

! possible standard names of precip from CMIP5, GRIB ECMWF and VECTRI dumps
  CHARACTER (LEN=6)  :: rainvar(4)=["precip","tp    ","pr    ","rain  "]  
  CHARACTER (LEN=11) :: tempvar(3) =["temperature","tmp        ","t2m        "]
  CHARACTER (LEN=31) :: time_units ! e.g. "days since 1990-11-25 00:00 UTC"

  CHARACTER (LEN = *), PARAMETER :: lat_name = "latitude"
  CHARACTER (LEN = *), PARAMETER :: lon_name = "longitude"
  CHARACTER (LEN = *), PARAMETER :: time_name = "time"
  CHARACTER (LEN = *), PARAMETER :: units = "units"
  CHARACTER (LEN = *), PARAMETER :: lat_units = "degrees_north"
  CHARACTER (LEN = *), PARAMETER :: lon_units = "degrees_east"

! logical for gridded runs
  LOGICAL :: l2d=.false.

! control of output file - default values.
  LOGICAL :: loutput_population=.true.  ! population
  LOGICAL :: loutput_rain=.true.        ! precipitation
  LOGICAL :: loutput_t2m=.true.         ! 2 metre temperature
  LOGICAL :: loutput_waterfrac=.true.   ! water fraction
  LOGICAL :: loutput_vector=.true.      ! vector density
  LOGICAL :: loutput_larvae=.true.      ! larvae density 
  LOGICAL :: loutput_lbiomass=.true.    ! larvae biomass
  LOGICAL :: loutput_pr=.false.         ! parasite ratio
  LOGICAL :: loutput_prd=.true.         ! detectable parasite ratio
  LOGICAL :: loutput_hbr=.true.         ! human bite rate
  LOGICAL :: loutput_cspr=.true.        ! circumsporozoite protein rate
  LOGICAL :: loutput_eir=.true.         ! Entomological Inoculation Rate
  LOGICAL :: loutput_immunity=.true.    ! Proportion of hosts immune
  LOGICAL :: loutput_cases=.true.       ! Proportion of new cases
  LOGICAL :: loutput_vecthostratio=.false. ! ratio of vectors to humans.

END MODULE mo_control
