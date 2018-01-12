MODULE mo_vectri
!-------------------------------------------------------- 
! VECTRI: VECtor borne disease infection model of TRIeste.
!
! Tompkins AM. 2011, ICTP
! tompkins@ictp.it
!
! main arrays and variables 
! for mosquitoes, larva -- all gridded data 
!--------------------------------------------------------
  USE mo_constants
  USE mo_control

  IMPLICIT NONE

  ! spatial dimensions
  INTEGER  :: nlat,nlon
 
  ! time dimensions 
  INTEGER :: nday        ! total integration time in days
  INTEGER :: start_date  ! start date yyyymmdd

  ! conversion of immune rate into tau parameter
  ! assumes 95% are immune at this value
  REAL :: rimmune_gain_tau

  ! daily death rate
  REAL :: rpop_death_rate_daily

  ! fraction of grid cell covered in water that offers breeding sites 
  ! this indicates pools, but not deep lakes
  ! There will be two components, permanent water bodies ponds
  ! and edges of lakes and rivers
  ! A second component that grows temporary pools with a 
  ! simple pond model, stored in rpuddle - 
  ! the total is in rwaterfrac *which is the prognostic variable*
  ! from v1.3.4 the permanent fraction is a function of space ready for 
  ! the LANDSAT/Sentinel classification maps.
  REAL, ALLOCATABLE :: rwaterfrac(:,:), rpuddle(:,:), rwaterperm(:,:) ! nlat,nlon          

  ! zoophilic rate, in places with low population rate zoophilic rate is 
  ! higher - Kileen et al. (2001) quotes two papers that show high trend 
  ! in Gambia, while sensitivity lower in Tanzania - will also depend
  ! on method since knockdown in huts could exaggerate anthropophilicity
  ! NOTE: zoophilic actually is anthropophilic rate!!! 1.0-zoo
  ! rbitezoo is just to make things faster...
  REAL, ALLOCATABLE :: rzoophilic(:,:), rbitezoo(:,:)  ! nx,ny
           
  ! number of female mosquitoes in gridcell 
  ! ngono defines the resolution of the gonotrophic cycle
  ! new vectors start in index 0 of ngono until they get a blood meal 
  !   they then move to box i, and advance according to degree days.
  ! likewise ninfv desribes the evolution of the parasite within the vector
  !   again, new vectors start and remain in box iinfv=zero until they 
  !   are infected during feeding.
  REAL, ALLOCATABLE :: rvect(:,:,:,:) ! rvect(0:ngono,0:ninfv,nx,ny)   

  ! larvae density per m**2 
  ! multiplied by water coverage to give total per box
  ! water coverage includes permanent pools, and adjusts with rains
  ! thus if 10 day rain rate is zero still have breeding sites
  REAL, ALLOCATABLE :: rlarv(:,:,:) ! rlarv(0:nlarv,nx,ny)     

  ! host infection state - will be replaced by nhost-space garki model
  REAL, ALLOCATABLE :: rhost(:,:,:,:) ! rhost(0:ninfh,nhost,nx,ny)   

! meteorological variables
  INTEGER, ALLOCATABLE  :: ndate(:)      ! ndate(nstep)         !          integer date
  REAL, ALLOCATABLE :: rtemp (:,:) ! rtemp(nstep,nx,ny)   ! (K)      T2m as a function of time and space
!  REAL, ALLOCATABLE :: rrh(:,:)    ! rrh(nstep,nx,ny)     ! (frac)   RH as a function of time and space
  REAL, ALLOCATABLE :: rrain(:,:)  ! rrain(nstep,nx,ny)   ! (mm/day) precip as a function of time and space
  REAL :: rrain_fillvalue, rtemp_fillvalue ! missing values in datasets

!----------------------------
! Non-gridded malaria arrays
!----------------------------
  REAL :: rmasslarv(0:nlarv)     ! mass of larvae 
  REAL :: rpdfvect2host(0:nbitepdf) ! PDF of transmission probability 
                                    ! as a function of bite number
  REAL :: rbiteSwgt(0:nbitepdf)     ! Susceptibles bite distribution weighting
                                    ! as a function of bite number
  REAL :: rbiteEIwgt(0:nbitepdf)    ! EI bite distribution weighting
                                    ! as a function of bite number
  REAL :: rsumrbiteSwgt, rsumrbiteEIwgt ! sum of the weightings

! rate constants for larvae maturation, array for convenience.
  REAL, DIMENSION (4,2) :: rlarvmature

! derived variable
  INTEGER ::  nrun ! length of run in timesteps


!--------------------------------
! control variables and constants
!--------------------------------
!
! climate data input type
!
! 1: constant input from namelist
! 2: station data from an ascii or netcdf file
! 3: gridded model/reanalysis/satellite data
!
  INTEGER  :: iclim_input

! station data file if required
  CHARACTER(len=200) :: station_file

! climate data file if required for gridded run
  CHARACTER(len=200) :: rainfile
  CHARACTER(len=200) :: tempfile

! population density, and population related parameters
  REAL :: rpopdensity2010, rpopdensity_FillValue
  REAL, ALLOCATABLE :: rpopdensity(:,:) ! rpopdensity(nx,ny)

! water occupancy rates for LUC experiments
  REAL, ALLOCATABLE :: rwateroccupancy(:,:) ! rwateroccupancy(nx,ny)

! region data
  REAL :: lon,dlon,lat,dlat
  REAL, ALLOCATABLE :: lons(:), lats(:) ! lon lat values

! --------------------------
! output diagnostics control
! --------------------------
  INTEGER :: ndaydiag ! diagnostics every n days
  INTEGER :: ncidout  ! ncdf file id for output
  INTEGER :: LonDimID, LatDimID, timeDimId

! ncids and indices for 2d output fields
  INTEGER, PARAMETER :: ncvar_ps=2
  INTEGER :: ncvar_rain(ncvar_ps)       
  INTEGER :: ncvar_t2m(ncvar_ps)       
  INTEGER :: ncvar_vector(ncvar_ps)    
  INTEGER :: ncvar_larvae(ncvar_ps)    
  INTEGER :: ncvar_lbiomass(ncvar_ps)   
  INTEGER :: ncvar_pr(ncvar_ps) 
!  INTEGER :: ncvar_vectinfect(ncvar_ps) 
  INTEGER :: ncvar_prd(ncvar_ps) 
!  INTEGER :: ncvar_cases(ncvar_ps)  
  INTEGER :: ncvar_eir(ncvar_ps)       
  INTEGER :: ncvar_cases(ncvar_ps)       
  INTEGER :: ncvar_immunity(ncvar_ps)       
  INTEGER :: ncvar_cspr(ncvar_ps)       
  INTEGER :: ncvar_hbr(ncvar_ps)      
  INTEGER :: ncvar_waterfrac(ncvar_ps)  
  INTEGER :: ncvar_vecthostratio(ncvar_ps)       
  INTEGER :: ncvar_popdensity(ncvar_ps)       
  INTEGER :: ndiag2d=0 ! number of diagnostics (defined in setup)
  REAL, ALLOCATABLE :: rdiag2d(:,:,:) ! nlon, nlat, ndiag2d 


! IDs for output variables
!  CHARACTER (len = *), PARAMETER :: temp_units = "celsius"
!  CHARACTER (len = *), PARAMETER :: rain_units = "mm/day"


END MODULE mo_vectri
